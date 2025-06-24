"""
Stages that implement GATK-gCNV.
"""

from functools import cache

import loguru
from cpg_flow import stage, targets, workflow
from cpg_utils import Path, config, hail_batch, to_path

from cpg_gcnv.jobs.AnnotateCnvsWithStrvctvre import annotate_cnvs_with_strvctvre
from cpg_gcnv.jobs.AnnotateCnvsWithSvAnnotate import queue_annotate_sv_jobs
from cpg_gcnv.jobs.AnnotateCohortCnv import submit_annotate_cohort_job
from cpg_gcnv.jobs.AnnotateDatasetCnv import submit_annotate_dataset_job
from cpg_gcnv.jobs.CallGermlineCnvsWithGatk import shard_gcnv
from cpg_gcnv.jobs.CollectReadCounts import collect_read_counts
from cpg_gcnv.jobs.DeterminePloidy import filter_and_determine_ploidy
from cpg_gcnv.jobs.FastCombineGCnvs import fast_merge_calls
from cpg_gcnv.jobs.JointSegmentCnvVcfs import run_joint_segmentation
from cpg_gcnv.jobs.MtToEsCnv import submit_es_job_for_dataset
from cpg_gcnv.jobs.PrepareIntervals import prepare_intervals
from cpg_gcnv.jobs.ProcessCohortCnvCallsToSgVcf import postprocess_unclustered_calls
from cpg_gcnv.jobs.RecalculateClusteredQuality import recalculate_clustered_calls
from cpg_gcnv.jobs.SplitAnnotatedVcfByDataset import split_mc_vcf_by_dataset
from cpg_gcnv.jobs.TrimOffSexChromosomes import trim_sex_chromosomes
from cpg_gcnv.jobs.UpgradePedWithInferredSex import upgrade_ped_file
from cpg_gcnv.utils import shard_items


@cache
def get_cohort_for_sgid(sgid: str) -> targets.Cohort:
    """
    Return the cohort that contains this sgid
    until a better central method exists this needs to run multiple times
    so build a method here with caching
    """
    for cohort_id in workflow.get_multicohort().get_cohorts():
        if sgid in cohort_id.get_sequencing_group_ids():
            return cohort_id
    raise ValueError(f'Could not find cohort for {sgid}')


@cache
def fixed_sg_order(cohort: targets.Cohort) -> list[str]:
    """
    using a method instead of a Stage - get the sorted list of IDs in this cohort

    Args:
        cohort ():

    Returns:

    """
    return sorted(cohort.get_sequencing_group_ids())


@stage.stage
class PrepareIntervals(stage.MultiCohortStage):
    """
    Interval preparation steps that don't require the sample read counts:
    PreprocessIntervals and AnnotateIntervals.
    This is a multicohort stage - we only ever co-process SGIDs on a unified capture
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path | str]:
        return {
            'preprocessed': self.prefix / 'preprocessed.interval_list',
            'annotated': self.prefix / 'annotated_intervals.tsv',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(multicohort)
        jobs = prepare_intervals(self.get_job_attrs(multicohort), outputs)
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=PrepareIntervals)
class CollectReadCounts(stage.SequencingGroupStage):
    """
    Per-sample stage that runs CollectReadCounts to produce .counts.tsv.gz files.
    """

    def expected_outputs(self, seqgroup: targets.SequencingGroup) -> dict[str, Path]:
        return {
            'counts': seqgroup.dataset.prefix() / 'gcnv' / f'{seqgroup.id}.counts.tsv.gz',
            'index': seqgroup.dataset.prefix() / 'gcnv' / f'{seqgroup.id}.counts.tsv.gz.tbi',
        }

    def queue_jobs(self, seqgroup: stage.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(seqgroup)

        if seqgroup.cram is None:
            raise ValueError(f'No CRAM file found for {seqgroup}')

        job = collect_read_counts(
            intervals_path=inputs.as_path(workflow.get_multicohort(), PrepareIntervals, 'preprocessed'),
            cram_path=seqgroup.cram,
            job_attrs=self.get_job_attrs(seqgroup),
            output_base_path=str(outputs['counts']).removesuffix('.counts.tsv.gz'),
        )
        return self.make_outputs(seqgroup, data=outputs, jobs=job)


@stage.stage(required_stages=[PrepareIntervals, CollectReadCounts])
class DeterminePloidy(stage.CohortStage):
    """
    The non-sharded cohort-wide gCNV steps after read counts have been collected:
    FilterIntervals and DetermineGermlineContigPloidy. These outputs represent
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        cohort_prefix = self.get_stage_cohort_prefix(cohort)
        return {
            'filtered': cohort_prefix / 'filtered.interval_list',
            'calls': cohort_prefix / 'ploidy-calls.tar.gz',
            'model': cohort_prefix / 'ploidy-model.tar.gz',
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(cohort)

        prep_intervals = inputs.as_dict(workflow.get_multicohort(), PrepareIntervals)

        # pull all per-sgid files from previous stage
        random_read_counts = inputs.as_path_by_target(CollectReadCounts, 'counts')

        # order those WRT the set ordering
        ordered_read_counts = [random_read_counts[seqgroup] for seqgroup in fixed_sg_order(cohort)]

        job = filter_and_determine_ploidy(
            ploidy_priors_path=config.config_retrieve(['references', 'gatk_sv', 'contig_ploidy_priors']),
            preprocessed_intervals_path=prep_intervals['preprocessed'],
            annotated_intervals_path=prep_intervals['annotated'],
            counts_paths=ordered_read_counts,
            job_attrs=self.get_job_attrs(cohort),
            output_paths=outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=job)


@stage.stage(required_stages=DeterminePloidy)
class UpgradePedWithInferred(stage.CohortStage):
    """
    Don't trust the metamist pedigrees, update with inferred sexes

    There's some hijinx to be done here - on the first pass this will never detect aneuploidies, because the files
    don't exist yet.

    We can iteratively run this pipeline - if aneuploid samples are detected in the later stages, they will be
    automatically detected by this step, or can be manually added from config
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        cohort_prefix = self.get_stage_cohort_prefix(cohort)
        return {
            'aneuploidy_samples': cohort_prefix / 'aneuploidies.txt',
            'pedigree': cohort_prefix / 'inferred_sex_pedigree.ped',
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(cohort)
        ploidy_inputs = hail_batch.get_batch().read_input(inputs.as_dict(cohort, DeterminePloidy)['calls'])
        ped_path = cohort.write_ped_file(self.get_stage_cohort_prefix(cohort, category='tmp') / 'pedigree.ped')
        job = upgrade_ped_file(
            ped_file=ped_path,
            new_output=outputs['pedigree'],
            aneuploidies=outputs['aneuploidy_samples'],
            ploidy_tar=ploidy_inputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=job)


@stage.stage(required_stages=[PrepareIntervals, CollectReadCounts, DeterminePloidy])
class CallGermlineCnvsWithGatk(stage.CohortStage):
    """
    The cohort-wide GermlineCNVCaller step, sharded across genome regions.
    This is separate from the DeterminePloidy stage so that the ProcessCohortCnvCallsToSgVcf
    stage can pick out this stage's sharded inputs easily.
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        return {name: self.get_stage_cohort_prefix(cohort) / f'{name}.tar.gz' for name in shard_items(name_only=True)}

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(cohort)
        determine_ploidy = inputs.as_dict(cohort, DeterminePloidy)
        prep_intervals = inputs.as_dict(workflow.get_multicohort(), PrepareIntervals)

        # pull all per-sgid files from previous stage
        random_read_counts = inputs.as_path_by_target(CollectReadCounts, 'counts')

        # order per-sgid files WRT the set ordering
        ordered_read_counts = [random_read_counts[seqgroup] for seqgroup in fixed_sg_order(cohort)]

        jobs = shard_gcnv(
            annotated_intervals_path=prep_intervals['annotated'],
            filtered_intervals_path=determine_ploidy['filtered'],
            ploidy_calls_path=determine_ploidy['calls'],
            counts_paths=ordered_read_counts,
            job_attrs=self.get_job_attrs(cohort),
            output_paths=outputs,
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=[DeterminePloidy, CallGermlineCnvsWithGatk])
class ProcessCohortCnvCallsToSgVcf(stage.SequencingGroupStage):
    """
    Produces final individual VCF results by running PostprocessGermlineCNVCalls.
    """

    def expected_outputs(self, seqgroup: targets.SequencingGroup) -> dict[str, Path]:
        """
        output paths here are per-SGID, but stored in the directory structure indicating the whole MCohort
        """

        # identify the cohort that contains this SGID
        this_cohort: targets.Cohort = get_cohort_for_sgid(seqgroup.id)

        # this job runs per sample, on results with a cohort context
        # so we need to write the outputs to a cohort-specific location
        cohort_prefix = self.get_stage_cohort_prefix(this_cohort)
        return {
            'intervals': cohort_prefix / f'{seqgroup.id}.intervals.vcf.gz',
            'intervals_index': cohort_prefix / f'{seqgroup.id}.intervals.vcf.gz.tbi',
            'segments': cohort_prefix / f'{seqgroup.id}.segments.vcf.gz',
            'segments_index': cohort_prefix / f'{seqgroup.id}.segments.vcf.gz.tbi',
            'ratios': cohort_prefix / f'{seqgroup.id}.ratios.tsv',
        }

    def queue_jobs(self, seqgroup: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(seqgroup)

        # identify the cohort that contains this SGID
        this_cohort = get_cohort_for_sgid(seqgroup.id)

        determine_ploidy = inputs.as_dict(this_cohort, DeterminePloidy)

        # pull the sgid ordering for this cohort
        sgid_ordering = fixed_sg_order(this_cohort)

        jobs = postprocess_unclustered_calls(
            ploidy_calls_path=determine_ploidy['calls'],
            shard_paths=inputs.as_dict(this_cohort, CallGermlineCnvsWithGatk),
            sample_index=sgid_ordering.index(seqgroup.id),
            job_attrs=self.get_job_attrs(seqgroup),
            output_prefix=str(self.get_stage_cohort_prefix(this_cohort) / seqgroup.id),
        )
        return self.make_outputs(seqgroup, data=outputs, jobs=jobs)


@stage.stage(required_stages=[ProcessCohortCnvCallsToSgVcf, UpgradePedWithInferred])
class TrimOffSexChromosomes(stage.CohortStage):
    """
    Trim off sex chromosomes for gCNV VCFs where the SGID is detected to be Aneuploid
    The dependency chain here is a number of CohortStages, i.e. thetargets.MultiCohort as a whole
    isn't relevant to determining aneuploidy. As a result we're happy writing this to a targets.Cohort-specific path
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path | str]:
        cohort_prefix = self.get_stage_cohort_prefix(cohort)

        sgids_in_this_cohort = cohort.get_sequencing_group_ids()

        # returning an empty dictionary might cause the pipeline setup to break?
        # returning a string won't cause an existince check on this file, but will trick the pipeline into running?
        return_dict: dict[str, Path | str] = {
            'placeholder': str(cohort_prefix / 'placeholder.txt'),
        }

        # load up the file of aneuploidies - I don't think the pipeline supports passing an input directly here
        # so... I'm making a similar path and manually string-replacing it
        aneuploidy_file = str(cohort_prefix / 'aneuploidies.txt').replace(
            self.name,
            'UpgradePedWithInferred',
        )

        # optionally pick up aneuploid samples from the config
        aneuploid_samples: list[str] = config.config_retrieve(['gCNV', 'aneuploid_samples'], [])

        if (aneuploidy_path := to_path(aneuploidy_file)).exists():
            # read the identified aneuploidy samples file
            with aneuploidy_path.open() as handle:
                # iterate over the lines
                for line in handle:
                    # find the SGID
                    sgid = line.strip()

                    # could be an empty newline
                    if not sgid or sgid not in sgids_in_this_cohort:
                        continue

                    aneuploid_samples.append(sgid)

        for sgid in set(aneuploid_samples):
            if sgid not in sgids_in_this_cohort:
                continue
            return_dict[sgid] = cohort_prefix / f'{sgid}.segments.vcf.bgz'

        if len(return_dict) > 1:
            # if don't need the placeholder, remove it
            return_dict.pop('placeholder')

        return return_dict

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        For each of the SGIDs which are identified as aneuploid, create a version
        with the X & Y chromosomes trimmed off
        Plan to generate every file, so that the stage can be forced to re-run if needed
        """
        expected = self.expected_outputs(cohort)
        germline_calls = inputs.as_dict_by_target(ProcessCohortCnvCallsToSgVcf)

        cohort_segment_vcfs = {sgid: germline_calls[sgid]['segments'] for sgid in cohort.get_sequencing_group_ids()}
        jobs = trim_sex_chromosomes(
            sgid_to_output=expected,
            segment_vcfs=cohort_segment_vcfs,
            job_attrs=self.get_job_attrs(cohort),
        )
        return self.make_outputs(cohort, data=expected, jobs=jobs)


@stage.stage(
    required_stages=[
        TrimOffSexChromosomes,
        ProcessCohortCnvCallsToSgVcf,
        PrepareIntervals,
        UpgradePedWithInferred,
    ],
)
class JointSegmentCnvVcfs(stage.CohortStage):
    """
    various config elements scavenged from https://github.com/broadinstitute/gatk/blob/cfd4d87ec29ac45a68f13a37f30101f326546b7d/scripts/cnv_cromwell_tests/germline/cnv_germline_case_scattered_workflow.json#L26
    continuing adaptation of https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/germline/joint_call_exome_cnvs.wdl
    takes the individual VCFs and runs the joint segmentation step
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path]:
        cohort_prefix = self.get_stage_cohort_prefix(cohort)
        return {
            'clustered_vcf': cohort_prefix / 'JointClusteredSegments.vcf.gz',
            'clustered_vcf_idx': cohort_prefix / 'JointClusteredSegments.vcf.gz.tbi',
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        So, this is a tricksy lil monster -
        Conducts a semi-heirarchical merge of the individual VCFs
        - First merge the segment files in blocks, to produce intermediate merges
        - Then merge those intermediate merges to produce the final result
        """

        outputs = self.expected_outputs(cohort)

        # get the individual Segment VCFs
        cnv_vcfs = inputs.as_dict_by_target(ProcessCohortCnvCallsToSgVcf)

        # and the dict of trimmed VCFs (can be empty)
        trimmed_vcfs = inputs.as_dict(cohort, TrimOffSexChromosomes)

        # pull the json file with the sgid ordering
        sgid_ordering = fixed_sg_order(cohort)

        # for each SGID, either get the sex chrom-trimmed one, or the default
        all_vcfs: list[str] = []
        for sgid in sgid_ordering:
            if sgid in trimmed_vcfs:
                loguru.logger.info(f'Using XY-trimmed VCF for {sgid}')
                all_vcfs.append(str(trimmed_vcfs[sgid]))
            elif sgid in cnv_vcfs:
                all_vcfs.append(str(cnv_vcfs[sgid]['segments']))
            else:
                raise ValueError(f'No VCF found for {sgid}')

        # get the intervals
        intervals = inputs.as_str(workflow.get_multicohort(), PrepareIntervals, 'preprocessed')

        pedigree = inputs.as_str(cohort, UpgradePedWithInferred, 'pedigree')

        jobs = run_joint_segmentation(
            segment_vcfs=all_vcfs,
            pedigree=pedigree,
            intervals=intervals,
            tmp_prefix=self.get_stage_cohort_prefix(cohort, category='tmp') / 'intermediate_jointseg',
            output_path=outputs['clustered_vcf'],
            job_attrs=self.get_job_attrs(cohort),
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(
    required_stages=[JointSegmentCnvVcfs, CallGermlineCnvsWithGatk, ProcessCohortCnvCallsToSgVcf, DeterminePloidy],
)
class RecalculateClusteredQuality(stage.SequencingGroupStage):
    """
    following joint segmentation, we need to post-process the clustered breakpoints
    this recalculates each sample's quality scores based on new breakpoints, and
    filters low QS or high AF calls
    https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/germline/joint_call_exome_cnvs.wdl#L113

    This is done as another pass through PostprocessGermlineCNVCalls, with prior/clustered results
    """

    def expected_outputs(self, seqgroup: targets.SequencingGroup) -> dict[str, Path]:
        # identify the cohort that contains this SGID
        this_cohort = get_cohort_for_sgid(seqgroup.id)

        cohort_prefix = self.get_stage_cohort_prefix(this_cohort)

        # this job runs per sample, on results with a cohort context
        # so we need to write the outputs to a cohort-specific location
        return {
            'genotyped_intervals_vcf': cohort_prefix / f'{seqgroup.id}.intervals.vcf.gz',
            'genotyped_intervals_vcf_index': cohort_prefix / f'{seqgroup.id}.intervals.vcf.gz.tbi',
            'genotyped_segments_vcf': cohort_prefix / f'{seqgroup.id}.segments.vcf.gz',
            'genotyped_segments_vcf_index': cohort_prefix / f'{seqgroup.id}.segments.vcf.gz.tbi',
            'denoised_copy_ratios': cohort_prefix / f'{seqgroup.id}.ratios.tsv',
            'qc_status_file': cohort_prefix / f'{seqgroup.id}.qc_status.txt',
        }

    def queue_jobs(self, sequencing_group: targets.SequencingGroup, inputs: stage.StageInput) -> stage.StageOutput:
        expected_out = self.expected_outputs(sequencing_group)

        # identify the cohort that contains this SGID
        this_cohort = get_cohort_for_sgid(sequencing_group.id)

        # get the clustered VCF from the previous stage
        joint_seg = inputs.as_dict(this_cohort, JointSegmentCnvVcfs)

        determine_ploidy = inputs.as_dict(this_cohort, DeterminePloidy)
        gcnv_call_inputs = inputs.as_dict(sequencing_group, ProcessCohortCnvCallsToSgVcf)

        # pull the json file with the sgid ordering
        sgid_ordering = fixed_sg_order(this_cohort)

        jobs = recalculate_clustered_calls(
            ploidy_calls_path=determine_ploidy['calls'],
            shard_paths=inputs.as_dict(this_cohort, CallGermlineCnvsWithGatk),
            sample_index=sgid_ordering.index(sequencing_group.id),
            job_attrs=self.get_job_attrs(sequencing_group),
            output_prefix=str(self.get_stage_cohort_prefix(this_cohort) / sequencing_group.id),
            clustered_vcf=str(joint_seg['clustered_vcf']),
            intervals_vcf=str(gcnv_call_inputs['intervals']),
            qc_file=str(expected_out['qc_status_file']),
        )
        return self.make_outputs(sequencing_group, data=expected_out, jobs=jobs)


@stage.stage(required_stages=RecalculateClusteredQuality)
class FastCombineGCNVs(stage.CohortStage):
    """
    Produces final multi-sample VCF results by running a merge
    """

    def expected_outputs(self, cohort: targets.Cohort) -> dict[str, Path | str]:
        """
        This is now explicitly continuing from multicohort work, so the output path must include
        pointers to both thetargets.MultiCohort and the targets.Cohort
        """
        cohort_prefix = self.get_stage_cohort_prefix(cohort)
        return {
            'combined_calls': cohort_prefix / 'gcnv_joint_call.vcf.bgz',
            'combined_calls_index': cohort_prefix / 'gcnv_joint_call.vcf.bgz.tbi',
        }

    def queue_jobs(self, cohort: targets.Cohort, inputs: stage.StageInput) -> stage.StageOutput:
        outputs = self.expected_outputs(cohort)
        gcnv_vcfs = inputs.as_dict_by_target(RecalculateClusteredQuality)
        all_vcfs = [str(gcnv_vcfs[sgid]['genotyped_segments_vcf']) for sgid in cohort.get_sequencing_group_ids()]

        jobs = fast_merge_calls(
            sg_vcfs=all_vcfs,
            job_attrs=self.get_job_attrs(cohort),
            output_path=outputs['combined_calls'],
        )
        return self.make_outputs(cohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=FastCombineGCNVs, analysis_type='cnv')
class MergeCohortsgCNV(stage.MultiCohortStage):
    """
    Takes all the per-targets.Cohort results and merges them into a pseudocallset
    We could use Jasmine for a better merge
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.prefix / 'multi_cohort_gcnv.vcf.bgz'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        re-uses the FastCombineGCNVs job to merge the targets.Cohort-level VCFs into atargets.MultiCohort VCF

        Args:
            multicohort ():
            inputs (stage.StageInput): link to FastCombineGCNVs outputs

        Returns:
            the bcftools merge job to join the targets.Cohort-VCFs into atargets.MultiCohort VCF
        """
        output = self.expected_outputs(multicohort)
        cohort_merges = inputs.as_dict_by_target(FastCombineGCNVs)
        cohort_vcfs = [str(cohort_merges[cohort.id]['combined_calls']) for cohort in multicohort.get_cohorts()]

        job_or_none = fast_merge_calls(
            sg_vcfs=cohort_vcfs,
            job_attrs=self.get_job_attrs(multicohort),
            output_path=output,
        )

        return self.make_outputs(multicohort, data=output, jobs=job_or_none)


@stage.stage(required_stages=MergeCohortsgCNV, analysis_type='cnv', analysis_keys=['annotated_vcf'])
class AnnotateCnvsWithSvAnnotate(stage.MultiCohortStage):
    """
    Smaller, direct annotation using SvAnnotate
    Add annotations, such as the inferred function and allele frequencies of variants, to final VCF.

    This is a full clone of the GATK-SV pipeline Cromwell stage, but use on a slightly
    different output. Trying to work out the best way to handle this through inheritance

    Annotations methods include:
    * Functional annotation - annotate SVs with inferred functional consequence on
      protein-coding regions, regulatory regions such as UTR and promoters, and other
      non-coding elements.
    * Allele frequency annotation - annotate SVs with their allele frequencies across
      all samples, and samples of specific sex, as well as specific subpopulations.
    * Allele Frequency annotation with external callset - annotate SVs with the allele
      frequencies of their overlapping SVs in another callset, e.g. gnomad SV callset.
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> dict[str, Path]:
        # kinda important to keep this as a dictionary for the pipeline - extracted by name as cromwell outputs
        return {
            'annotated_vcf': self.prefix / 'merged_gcnv_annotated.vcf.bgz',
            'annotated_vcf_index': self.prefix / 'merged_gcnv_annotated.vcf.bgz.tbi',
        }

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        configure and queue jobs for SV annotation
        passing the VCF Index has become implicit, which may be a problem for us
        """
        outputs = self.expected_outputs(multicohort)

        input_vcf = inputs.as_str(multicohort, MergeCohortsgCNV)

        jobs = queue_annotate_sv_jobs(
            multicohort=multicohort,
            prefix=self.prefix,
            input_vcf=input_vcf,
            outputs=outputs,
        )
        return self.make_outputs(multicohort, data=outputs, jobs=jobs)


@stage.stage(required_stages=AnnotateCnvsWithSvAnnotate, analysis_type='cnv')
class AnnotateCnvsWithStrvctvre(stage.MultiCohortStage):
    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.prefix / 'cnv_strvctvre_annotated.vcf.bgz'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(multicohort)

        input_vcf = inputs.as_str(multicohort, AnnotateCnvsWithSvAnnotate, key='annotated_vcf')

        job = annotate_cnvs_with_strvctvre(
            input_vcf=input_vcf,
            output_vcf=str(output),
            job_attrs=self.get_job_attrs() | {'tool': 'strvctvre'},
        )

        return self.make_outputs(multicohort, data=output, jobs=job)


@stage.stage(required_stages=AnnotateCnvsWithStrvctvre, analysis_type='single_dataset_cnv_annotated')
class SplitAnnotatedCnvVcfByDataset(stage.DatasetStage):
    """
    takes the whole MultiCohort annotated VCF
    splits it up into separate VCFs for each dataset
    """

    def expected_outputs(self, dataset: targets.dataset) -> Path:
        return (
            dataset.prefix()
            / workflow.get_workflow().name
            / workflow.get_workflow().output_version
            / self.name
            / 'annotated_cnv.vcf.bgz'
        )

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        output = self.expected_outputs(dataset)
        input_vcf = hail_batch.get_batch().read_input(
            inputs.as_path(workflow.get_multicohort(), AnnotateCnvsWithStrvctvre),
        )

        job = split_mc_vcf_by_dataset(
            input_vcf=input_vcf,
            dataset=dataset,
            output=str(output),
            job_attrs=self.get_job_attrs(dataset),
        )

        return self.make_outputs(dataset, data=output, jobs=job)


@stage.stage(required_stages=AnnotateCnvsWithStrvctvre, analysis_type='cnv')
class AnnotateCohortCnv(stage.MultiCohortStage):
    """
    Rearrange the annotations across the cohort to suit Seqr
    """

    def expected_outputs(self, multicohort: targets.MultiCohort) -> Path:
        return self.prefix / 'gcnv_annotated_cohort.mt'

    def queue_jobs(self, multicohort: targets.MultiCohort, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Fire up the job to ingest the cohort VCF as a MT, and rearrange the annotations
        """

        vcf_path = inputs.as_str(target=multicohort, stage=AnnotateCnvsWithStrvctvre)

        output = self.expected_outputs(multicohort)

        job = submit_annotate_cohort_job(
            input_vcf=vcf_path,
            output_mt=output,
            checkpoint=self.tmp_prefix / 'checkpoint.mt',
            attributes=self.get_job_attrs(multicohort),
        )

        return self.make_outputs(multicohort, data=output, jobs=job)


@stage.stage(required_stages=AnnotateCohortCnv, analysis_type='cnv')
class AnnotateDatasetCnv(stage.DatasetStage):
    """
    Subset the MT to be this targets.Dataset only, then work up all the genotype values
    """

    def expected_outputs(self, dataset: targets.Dataset) -> Path:
        """
        Expected to generate a matrix table
        """
        return dataset.prefix() / 'mt' / f'gCNV-{workflow.get_workflow().output_version}-{dataset.name}.mt'

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Subsets the multicohort MT to this dataset only, then brings genotype data into row annotations

        Args:
            dataset ():
            inputs ():
        """

        mt_in = inputs.as_str(target=workflow.get_multicohort(), stage=AnnotateCohortCnv)

        output = self.expected_outputs(dataset)

        job = submit_annotate_dataset_job(input_mt=mt_in, output_mt=output, attributes=self.get_job_attrs(dataset))

        return self.make_outputs(dataset, data=output, jobs=job)


@stage.stage(
    required_stages=[AnnotateDatasetCnv],
    analysis_type='es-index',
    analysis_keys=['index_name'],
    # https://github.com/populationgenomics/metamist/issues/539
    update_analysis_meta=lambda x: {'seqr-dataset-type': 'CNV'},  # noqa: ARG005
)
class MtToEsCnv(stage.DatasetStage):
    """
    Create a Seqr index
    """

    def expected_outputs(self, dataset: targets.Dataset) -> dict[str, str | Path]:
        """
        Expected to generate a Seqr index, which is not a file
        """
        index_name = f'{dataset.name}-exome-gCNV-{workflow.get_workflow().run_timestamp}'.lower()
        return {
            'index_name': index_name,
            'done_flag': dataset.prefix() / 'es' / f'{index_name}.done',
        }

    def queue_jobs(self, dataset: targets.Dataset, inputs: stage.StageInput) -> stage.StageOutput:
        """
        Freshly liberated from the clutches of DataProc
        Uses the script cpg_workflows/dataproc_scripts/mt_to_es_free_of_dataproc.py
        The script was registered in setup.py with a console entrypoint
        This requires a code version >= 1.25.14 in the worker job image to operate

        gCNV indexes have never been spotted over a GB in the wild, so we use minimum storage

        Args:
            dataset (targets.Dataset):
            inputs ():
        """

        # get the absolute path to the MT
        mt_path = inputs.as_str(target=dataset, stage=AnnotateDatasetCnv)

        outputs = self.expected_outputs(dataset)

        job = submit_es_job_for_dataset(
            mt_path=mt_path,
            index_name=outputs['index_name'],
            done_flag=outputs['done_flag'],
            dataset=dataset.name,
        )

        return self.make_outputs(dataset, data=outputs, jobs=job)
