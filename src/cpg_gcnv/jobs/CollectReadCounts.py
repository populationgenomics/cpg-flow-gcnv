from typing import TYPE_CHECKING

from cpg_flow.resources import HIGHMEM
from cpg_utils.config import image_path
from cpg_utils.hail_batch import fasta_res_group, get_batch

if TYPE_CHECKING:
    from cpg_flow.filetypes import CramPath
    from cpg_utils import Path
    from hailtop.batch.job import Job


def collect_read_counts(
    intervals_path: 'Path',
    cram_path: 'CramPath',
    job_attrs: dict[str, str],
    output_base_path: str,
) -> 'Job':
    job = get_batch().new_job('Collect gCNV read counts', job_attrs | {'tool': 'gatk CollectReadCounts'})
    job.image(image_path('gatk_gcnv'))

    reference = fasta_res_group(get_batch())

    job.declare_resource_group(
        counts={
            'tsv.gz': '{root}.counts.tsv.gz',
            'tsv.gz.tbi': '{root}.counts.tsv.gz.tbi',
        },
    )

    # set highmem resources for this job
    job_res = HIGHMEM.request_resources(ncpu=2, storage_gb=10)
    job_res.set_to_job(job)

    job.cpu(2).storage('10Gi').memory('highmem')

    job.command(f"""
    gatk --java-options "{job_res.java_mem_options()}" \
        CollectReadCounts \
        --reference {reference.base} \
        --intervals {intervals_path} \
        --interval-merging-rule OVERLAPPING_ONLY \
        --input {cram_path.path} \
        --read-index {cram_path.index_path} \
        --format TSV \
        --output {job.counts}.counts.tsv
        """)

    job.command(f'bgzip {job.counts}.counts.tsv')
    job.command(f'gatk IndexFeatureFile --input {job.counts["counts.tsv.gz"]}')
    get_batch().write_output(job.counts, output_base_path)
    return job
