import re
from typing import TYPE_CHECKING

from cpg_flow.resources import HIGHMEM
from cpg_utils import to_path
from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import authenticate_cloud_credentials_in_job, fasta_res_group, get_batch

if TYPE_CHECKING:
    from cpg_flow.targets import MultiCohort
    from cpg_utils import Path
    from hailtop.batch.job import BashJob


PED_FAMILY_ID_REGEX = re.compile(r'(^[A-Za-z0-9_]+$)')


def chunks(iterable, chunk_size):
    """
    Yield successive n-sized chunks from an iterable

    Args:
        iterable (): any iterable - tuple, str, list, set
        chunk_size (): size of intervals to return

    Returns:
        intervals of requested size across the collection
    """

    if isinstance(iterable, set):
        iterable = list(iterable)

    for i in range(0, len(iterable), chunk_size):
        yield iterable[i : (i + chunk_size)]


def clean_ped_family_ids(ped_line: str) -> str:
    """
    Takes each line in the pedigree and cleans it up
    If the family ID already conforms to expectations, no action
    If the family ID fails, replace all non-alphanumeric/non-underscore
    characters with underscores

    >>> clean_ped_family_ids('family1\tchild1\t0\t0\t1\t0\\n')
    'family1\tchild1\t0\t0\t1\t0\\n'
    >>> clean_ped_family_ids('family-1-dirty\tchild1\t0\t0\t1\t0\\n')
    'family_1_dirty\tchild1\t0\t0\t1\t0\\n'

    Args:
        ped_line (str): line from the pedigree file, unsplit

    Returns:
        the same line with a transformed family id
    """

    split_line = ped_line.rstrip().split('\t')

    if re.match(PED_FAMILY_ID_REGEX, split_line[0]):
        return ped_line

    # if the family id is not valid, replace failing characters with underscores
    split_line[0] = re.sub(r'[^A-Za-z0-9_]', '_', split_line[0])

    # return the rebuilt string, with a newline at the end
    return '\t'.join(split_line) + '\n'


def make_combined_ped(cohort: 'MultiCohort', prefix: 'Path') -> 'Path':
    """
    Create cohort + ref panel PED.
    Concatenating all samples across all datasets with ref panel

    See #578 - there are restrictions on valid characters in PED file
    """
    combined_ped_path = prefix / 'ped_with_ref_panel.ped'
    conf_ped_path = config_retrieve(['references', 'broad', 'ped_file'])
    with combined_ped_path.open('w') as out:
        with cohort.write_ped_file().open() as f:
            # layer of family ID cleaning
            for line in f:
                out.write(clean_ped_family_ids(line))
        # The ref panel PED doesn't have any header, so can safely concatenate:
        with to_path(conf_ped_path).open() as f:
            out.write(f.read())
    return combined_ped_path


def postprocess_calls(
    ploidy_calls_path: 'Path',
    shard_paths: dict[str, 'Path'],
    sample_index: int,
    job_attrs: dict[str, str],
    output_prefix: str,
    clustered_vcf: str | None = None,
    intervals_vcf: str | None = None,
    qc_file: str | None = None,
) -> 'BashJob':
    if any([clustered_vcf, intervals_vcf, qc_file]) and not all([clustered_vcf, intervals_vcf, qc_file]):
        raise ValueError(
            'If any of clustered_vcf, intervals_vcf, or qc_file are provided, all must be provided',
        )

    job_name = 'Postprocess gCNV calls with clustered VCF' if clustered_vcf else 'Postprocess gCNV calls'

    job = get_batch().new_bash_job(job_name, job_attrs | {'tool': 'gatk PostprocessGermlineCNVCalls'})
    job.image(image_path('gatk_gcnv'))

    # set highmem resources for this job
    job_res = HIGHMEM.request_resources(ncpu=2, storage_gb=20)
    job_res.set_to_job(job)
    authenticate_cloud_credentials_in_job(job)

    reference = fasta_res_group(get_batch())

    ploidy_calls_tarball = get_batch().read_input(ploidy_calls_path)

    job.command(f'tar -xzf {ploidy_calls_tarball} -C $BATCH_TMPDIR/inputs')

    model_shard_args = ''
    calls_shard_args = ''

    for name, path in [(shard, shard_paths[shard]) for shard in shard_items(name_only=True)]:
        shard_tar = get_batch().read_input(path)
        job.command(f'tar -xzf {shard_tar} -C $BATCH_TMPDIR/inputs')
        model_shard_args += f' --model-shard-path $BATCH_TMPDIR/inputs/{name}-model'
        calls_shard_args += f' --calls-shard-path $BATCH_TMPDIR/inputs/{name}-calls'

    allosomal_contigs_args = ' '.join(
        [f'--allosomal-contig {c}' for c in config_retrieve(['workflow', 'allosomal_contigs'], [])]
    )

    # declare all output files in advance
    job.declare_resource_group(
        output={
            'intervals.vcf.gz': '{root}/intervals.vcf.gz',
            'intervals.vcf.gz.tbi': '{root}/intervals.vcf.gz.tbi',
            'segments.vcf.gz': '{root}/segments.vcf.gz',
            'segments.vcf.gz.tbi': '{root}/segments.vcf.gz.tbi',
            'ratios.tsv': '{root}/ratios.tsv',
        },
    )

    extra_args = ''
    if clustered_vcf:
        local_clusters = get_batch().read_input_group(vcf=clustered_vcf, index=f'{clustered_vcf}.tbi').vcf
        local_intervals = get_batch().read_input_group(vcf=intervals_vcf, index=f'{intervals_vcf}.tbi').vcf
        extra_args += f"""--clustered-breakpoints {local_clusters} \
         --input-intervals-vcf {local_intervals} \
          -R {reference.base}
        """

    job.command(
        f"""
    OUTS=$(dirname {job.output['intervals.vcf.gz']})
    BATCH_OUTS=$(dirname $OUTS)
    mkdir $OUTS
    gatk --java-options "{job_res.java_mem_options()}" PostprocessGermlineCNVCalls \
      --sequence-dictionary {reference.dict} {allosomal_contigs_args} \
      --contig-ploidy-calls $BATCH_TMPDIR/inputs/ploidy-calls \
      {model_shard_args} {calls_shard_args} \
      --sample-index {sample_index} \
      --output-genotyped-intervals {job.output['intervals.vcf.gz']} \
      --output-genotyped-segments {job.output['segments.vcf.gz']} \
      --output-denoised-copy-ratios {job.output['ratios.tsv']} {extra_args}
    """,
    )

    # index the output VCFs
    job.command(f'tabix -f {job.output["intervals.vcf.gz"]}')
    job.command(f'tabix -f {job.output["segments.vcf.gz"]}')

    if clustered_vcf:
        max_events = config_retrieve(['workflow', 'gncv_max_events'])
        max_pass_events = config_retrieve(['workflow', 'gncv_max_pass_events'])

        # do some additional QC to determine pass/fail
        job.command(
            f"""
        #use awk instead of grep - grep returning no lines causes a pipefailure
        NUM_SEGMENTS=$(zcat {job.output['segments.vcf.gz']} | \
          awk '!/^#/ && !/0\\/0/ && !/\t0:1:/ {{count++}} END {{print count}}')
        NUM_PASS_SEGMENTS=$(zcat {job.output['segments.vcf.gz']} | \
          awk '!/^#/ && !/0\\/0/ && !/\t0:1:/ && /PASS/ {{count++}} END {{print count}}')
        if [ $NUM_SEGMENTS -lt {max_events} ]; then
            if [ $NUM_PASS_SEGMENTS -lt {max_pass_events} ]; then
              echo "PASS" >> {job.qc_file}
            else
              echo "EXCESSIVE_NUMBER_OF_PASS_EVENTS" >> {job.qc_file}
            fi
        else
            echo "EXCESSIVE_NUMBER_OF_EVENTS" >> {job.qc_file}
        fi
        cat {job.qc_file}
        """,
        )
        get_batch().write_output(job.qc_file, qc_file)

    get_batch().write_output(job.output, output_prefix)

    return job


def counts_input_args(counts_paths: list['Path']) -> str:
    args = ''
    for f in counts_paths:
        counts = get_batch().read_input_group(
            **{
                'counts.tsv.gz': str(f),
                'counts.tsv.gz.tbi': str(f) + '.tbi',
            },
        )['counts.tsv.gz']
        args += f' --input {counts}'

    return args


def shard_items(name_only: bool):
    if shard_partition := config_retrieve(['workflow', 'interval_shards']):
        count = len(shard_partition)
        for i, interval in enumerate(shard_partition, start=1):
            name = 'part' + 'c'.join([s.removeprefix('chr') for s in interval])
            if name_only:
                yield name
            else:
                chroms = '|'.join(interval)
                yield name, i, count, f"awk '/^@/ || $1 ~ /^({chroms})$/'"
    else:
        raise NotImplementedError('gCNV sharding technique not specified')
