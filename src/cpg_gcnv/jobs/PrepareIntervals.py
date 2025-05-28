from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def prepare_intervals(
    job_attrs: dict[str, str],
    output_paths: dict[str, Path],
) -> 'BashJob':
    """
    generate a job to prepare the intervals file for use by gCNV

    Args:
        job_attrs (dict): any params to attach to the job
        output_paths (dict): paths to write the output files to

    Returns:
        job (BashJob): the job object
    """
    job = hail_batch.get_batch().new_bash_job(
        'Prepare intervals',
        job_attrs | {'tool': 'gatk PreprocessIntervals/AnnotateIntervals'},
    )
    job.image(config.config_retrieve(['images', 'gatk_gcnv']))

    reference = hail_batch.fasta_res_group(hail_batch.get_batch())

    exclude_intervals_args = ' '.join(
        [
            f'--exclude-intervals {i}'
            for i in config.config_retrieve(
                ['workflow', 'exclude_intervals'],
                [],
            )
        ],
    )

    intervals = hail_batch.get_batch().read_input(config.config_retrieve(['workflow', 'intervals_path']))

    # give the file an accurate extension
    job.preprocessed.add_extension('.interval_list')

    job.command(f"""
    gatk PreprocessIntervals \\
        --reference {reference.base} \\
        --intervals {intervals} \\
        {exclude_intervals_args} \\
        --padding 250 \\
        --bin-length 0 \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --output {job.preprocessed}
    """)
    job.command(f"""
    gatk AnnotateIntervals \\
        --reference {reference.base} \\
        --intervals {job.preprocessed} \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --output {job.annotated}
    """)

    for key, path in output_paths.items():
        hail_batch.get_batch().write_output(job[key], path)

    return job
