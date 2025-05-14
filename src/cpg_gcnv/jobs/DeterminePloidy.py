from typing import TYPE_CHECKING

from cpg_flow.resources import HIGHMEM
from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch

if TYPE_CHECKING:
    from cpg_utils import Path
    from hailtop.batch.job import BashJob


def counts_input_getter(counts_paths: list['Path']) -> str:
    args = ''
    for f in counts_paths:
        counts = get_batch().read_input_group(
            **{
                'counts.tsv.gz': str(f),
                'counts.tsv.gz.tbi': str(f) + '.tbi',
            },
        )
        args += f' --input {counts["counts.tsv.gz"]}'

    return args


def filter_and_determine_ploidy(
    ploidy_priors_path: str,
    preprocessed_intervals_path: 'Path',
    annotated_intervals_path: 'Path',
    counts_paths: list['Path'],
    job_attrs: dict[str, str],
    output_paths: dict[str, 'Path'],
) -> 'BashJob':
    job = get_batch().new_bash_job(
        'Filter intervals and determine ploidy',
        job_attrs
        | {
            'tool': 'gatk FilterIntervals/DetermineGermlineContigPloidy',
        },
    )
    job.image(image_path('gatk_gcnv'))

    # set highmem resources for this job
    job_res = HIGHMEM.request_resources(ncpu=2, storage_gb=10)
    job_res.set_to_job(job)

    counts_input_args = counts_input_getter(counts_paths)

    preprocessed_intervals = get_batch().read_input(preprocessed_intervals_path)
    annotated_intervals = get_batch().read_input(annotated_intervals_path)

    job.filtered.add_extension('.interval_list')

    job.command(f"""
    gatk --java-options "{job_res.java_mem_options()}" FilterIntervals \
      --interval-merging-rule OVERLAPPING_ONLY \
      --intervals {preprocessed_intervals} --annotated-intervals {annotated_intervals} \
      {counts_input_args} \
      --output {job.filtered}
    """)

    # (Other arguments may be cloud URLs, but this *must* be a local file)
    ploidy_priors = get_batch().read_input(ploidy_priors_path)

    job.command(f"""
    gatk \
        --java-options "{job_res.java_mem_options()}" \
        DetermineGermlineContigPloidy \
        --interval-merging-rule OVERLAPPING_ONLY \
        --intervals {job.filtered} \
        --contig-ploidy-priors {ploidy_priors} \
        {counts_input_args} \
        --output $BATCH_TMPDIR \
        --output-prefix ploidy

    tar -czf {job.calls} -C $BATCH_TMPDIR ploidy-calls
    tar -czf {job.model} -C $BATCH_TMPDIR ploidy-model
    """)

    for key, path in output_paths.items():
        get_batch().write_output(job[key], path)
    return job
