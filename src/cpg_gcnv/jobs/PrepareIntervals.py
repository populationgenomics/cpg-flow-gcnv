from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import fasta_res_group, get_batch

if TYPE_CHECKING:
    from cpg_utils import Path
    from hailtop.batch.job import Job


def prepare_intervals(job_attrs: dict[str, str], output_paths: dict[str, 'Path']) -> 'Job':
    job = get_batch().new_job(
        'Prepare intervals',
        job_attrs | {'tool': 'gatk PreprocessIntervals/AnnotateIntervals'},
    )
    job.image(image_path('gatk_gcnv'))

    reference = fasta_res_group(get_batch())

    exclude_intervals_args = ' '.join(
        [f'--exclude-intervals {i}' for i in config_retrieve(['workflow', 'exclude_intervals'], [])]
    )

    intervals = get_batch().read_input(config_retrieve(['workflow', 'intervals_path']))
    job.command(f"""
    gatk PreprocessIntervals \
        --reference {reference.base} \
        --intervals {intervals} \
        {exclude_intervals_args} \
        --padding 250 \
        --bin-length 0 \
        --interval-merging-rule OVERLAPPING_ONLY \
        --output {job.preprocessed}
    """)

    # give the file an accurate extension
    job.preprocessed.add_extension('.interval_list')

    job.command(f"""
    gatk AnnotateIntervals \
        --reference {reference.base} \
        --intervals {job.preprocessed} \
        --interval-merging-rule OVERLAPPING_ONLY \
        --output {job.annotated}
    """)

    for key, path in output_paths.items():
        get_batch().write_output(job[key], str(path))
    return job
