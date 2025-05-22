from random import randint
from typing import TYPE_CHECKING

from cpg_flow import targets
from cpg_utils import Path, config, cromwell, hail_batch

from cpg_gcnv import utils

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


GATK_SV_COMMIT = 'dc145a52f76a6f425ac3f481171040e78c0cfeea'
ANNOTATION_WORKFLOW = 'AnnotateVcf'


def queue_annotate_sv_jobs(
    multicohort: targets.MultiCohort,
    prefix: Path,
    input_vcf: Path,
    outputs: dict[str, Path],
) -> 'list[BashJob]':
    """
    Create an Annotation job in the Cromwell workflow engine
    """

    input_dict: dict = {
        'vcf': input_vcf,
        'prefix': multicohort.name,
        'ped_file': utils.make_combined_ped(multicohort, prefix),
        'sv_per_shard': 5000,
        'external_af_population': config.config_retrieve(['references', 'gatk_sv', 'external_af_population']),
        'external_af_ref_prefix': config.config_retrieve(['references', 'gatk_sv', 'external_af_ref_bed_prefix']),
        'external_af_ref_bed': config.config_retrieve(['references', 'gnomad_sv']),
        'use_hail': False,
        'noncoding_bed': config.config_retrieve(['references', 'gatk_sv', 'noncoding_bed']),
        'protein_coding_gtf': config.config_retrieve(['references', 'gatk_sv', 'protein_coding_gtf']),
        'contig_list': config.config_retrieve(['references', 'broad', 'primary_contigs_list']),
        # images
        'gatk_docker': config.config_retrieve(['images', 'gatk_docker']),
        'sv_pipeline_docker': config.config_retrieve(['images', 'sv_pipeline_docker']),
        'sv_base_mini_docker': config.config_retrieve(['images', 'sv_base_mini_docker']),
    }

    # obtain upper and lower polling bounds for this job size
    polling_minimum = randint(40, 80)  # noqa: S311
    polling_maximum = randint(400, 800)  # noqa: S311

    # If a config section exists for this workflow, apply overrides
    if override := config.config_retrieve(['resource_overrides', 'AnnotateVcf'], False):
        input_dict |= override

    # Where Cromwell writes the output - Will be different from paths in expected_out_dict:
    output_prefix = f'gatk_sv/output/AnnotateVcf/{multicohort.analysis_dataset.name}'

    outputs_to_collect = {
        key: cromwell.CromwellOutputType.single_path(f'{ANNOTATION_WORKFLOW}.{key}') for key, value in outputs.items()
    }

    # pre-process input_dict
    paths_as_strings: dict = {}
    for key, value in input_dict.items():
        if isinstance(value, Path):
            paths_as_strings[f'{ANNOTATION_WORKFLOW}.{key}'] = str(value)
        elif isinstance(value, list | set):
            paths_as_strings[f'{ANNOTATION_WORKFLOW}.{key}'] = [str(v) for v in value]
        else:
            paths_as_strings[f'{ANNOTATION_WORKFLOW}.{key}'] = value

    submit_j, output_dict = cromwell.run_cromwell_workflow_from_repo_and_get_outputs(
        b=hail_batch.get_batch(),
        job_prefix=f'Annotate CNVs: {multicohort.analysis_dataset.name}',
        dataset=config.config_retrieve(['workflow', 'dataset']),
        repo='gatk-sv',
        commit=GATK_SV_COMMIT,
        cwd='wdl',
        workflow='AnnotateVcf.wdl',
        libs=['.'],
        output_prefix=output_prefix,
        input_dict=paths_as_strings,
        outputs_to_collect=outputs_to_collect,
        driver_image=config.config_retrieve(['workflow', 'driver_image']),
        copy_outputs_to_gcp=config.config_retrieve(['workflow', 'copy_outputs'], False),
        labels={'stage': 'AnnotateCnvsWithSvAnnotate', config.AR_GUID_NAME: config.try_get_ar_guid()},
        min_watch_poll_interval=polling_minimum,
        max_watch_poll_interval=polling_maximum,
    )

    copy_j = hail_batch.get_batch().new_job(f'Annotate CNVs: {multicohort.analysis_dataset.name}: copy outputs')
    hail_batch.authenticate_cloud_credentials_in_job(copy_j)
    copy_j.image(config.config_retrieve(['workflow', 'driver_image']))
    for key, resource in output_dict.items():
        copy_j.command(f'gcloud storage cp "$(cat {resource})" "{outputs[key]}"')
    return [submit_j, copy_j]
