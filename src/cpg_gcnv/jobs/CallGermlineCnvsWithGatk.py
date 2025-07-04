from typing import TYPE_CHECKING

from cpg_flow import resources
from cpg_flow import utils as cpg_utils
from cpg_utils import Path, config, hail_batch

from cpg_gcnv.utils import counts_input_args, shard_items

if TYPE_CHECKING:
    from hailtop.batch.job import Job


def shard_gcnv(
    annotated_intervals_path: Path,
    filtered_intervals_path: Path,
    ploidy_calls_path: Path,
    counts_paths: list[Path],
    job_attrs: dict[str, str],
    output_paths: dict[str, Path],
) -> 'list[Job]':
    batch_instance = hail_batch.get_batch()
    annotated_intervals = batch_instance.read_input(annotated_intervals_path)
    filtered_intervals = batch_instance.read_input(filtered_intervals_path)
    ploidy_calls_tarball = batch_instance.read_input(ploidy_calls_path)
    counts_input = counts_input_args(counts_paths)

    jobs: list[Job] = []

    for name, i, n, select_cmd in shard_items(name_only=False):
        if cpg_utils.can_reuse(output_paths[name]):
            continue

        job = batch_instance.new_job(
            f'Call germline CNVs shard {i} of {n}',
            job_attrs
            | {
                'tool': 'gatk GermlineCNVCaller',
            },
        )
        job.image(config.config_retrieve(['images', 'gatk_gcnv']))

        # set highmem resources for this job
        job_res = resources.HIGHMEM.request_resources(ncpu=8, mem_gb=52, storage_gb=10)
        job_res.set_to_job(job)

        job.command(f"""
        tar -xzf {ploidy_calls_tarball} -C $BATCH_TMPDIR

        {select_cmd} < {filtered_intervals} > job.shard.interval_list

        gatk \\
          --java-options "{job_res.java_mem_options()}" GermlineCNVCaller \\
          --run-mode COHORT --interval-merging-rule OVERLAPPING_ONLY \\
          --intervals job.shard.interval_list --annotated-intervals {annotated_intervals} \\
          {counts_input} \\
          --contig-ploidy-calls $BATCH_TMPDIR/ploidy-calls \\
          --output $BATCH_TMPDIR --output-prefix {name}

        tar -czf {job.shard_tarball} -C $BATCH_TMPDIR {name}-calls {name}-model {name}-tracking
        """)
        batch_instance.write_output(job.shard_tarball, output_paths[name])
        jobs.append(job)

    return jobs
