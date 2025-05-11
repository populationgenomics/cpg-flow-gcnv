from typing import TYPE_CHECKING
from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch
from cpg_flow.utils import can_reuse
from cpg_flow.resources import HIGHMEM

from cpg_gcnv.utils import counts_input_args, shard_items

if TYPE_CHECKING:
    from cpg_utils import Path
    from hailtop.batch.job import Job


def shard_gcnv(
    annotated_intervals_path: 'Path',
    filtered_intervals_path: 'Path',
    ploidy_calls_path: 'Path',
    counts_paths: 'list[Path]',
    job_attrs: dict[str, str],
    output_paths: 'dict[str, Path]',
) -> 'list[Job]':

    annotated_intervals = get_batch().read_input(annotated_intervals_path)
    filtered_intervals = get_batch().read_input(filtered_intervals_path)
    ploidy_calls_tarball = get_batch().read_input(ploidy_calls_path)
    counts_input = counts_input_args(counts_paths)

    jobs: list[Job] = []

    for name, i, n, select_cmd in shard_items(name_only=False):
        if can_reuse(output_paths[name]):
            continue

        job = get_batch().new_job(
            f'Call germline CNVs shard {i} of {n}',
            job_attrs
            | {
                'tool': 'gatk GermlineCNVCaller',
            },
        )
        job.image(image_path('gatk_gcnv'))

        # set highmem resources for this job
        job_res = HIGHMEM.request_resources(ncpu=8, mem_gb=52, storage_gb=10)
        job_res.set_to_job(job)

        job.command(f"""
        tar -xzf {ploidy_calls_tarball} -C $BATCH_TMPDIR

        {select_cmd} < {filtered_intervals} > {job.shard_intervals}

        gatk \
          --java-options "{job_res.java_mem_options()}" GermlineCNVCaller \
          --run-mode COHORT --interval-merging-rule OVERLAPPING_ONLY \
          --intervals {job.shard_intervals} --annotated-intervals {annotated_intervals} \
          {counts_input} \
          --contig-ploidy-calls $BATCH_TMPDIR/ploidy-calls \
          --output $BATCH_TMPDIR --output-prefix {name}

        tar -czf {job.shard_tarball} -C $BATCH_TMPDIR {name}-calls {name}-model {name}-tracking
        """)
        get_batch().write_output(job.shard_tarball, output_paths[name])
        jobs.append(job)

    return jobs
