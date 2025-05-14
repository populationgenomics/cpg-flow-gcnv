from typing import TYPE_CHECKING

from cpg_flow.resources import HIGHMEM
from cpg_utils.config import image_path
from cpg_utils.hail_batch import fasta_res_group, get_batch

if TYPE_CHECKING:
    from cpg_flow.filetypes import CramPath
    from cpg_utils import Path
    from hailtop.batch.job import BashJob


def collect_read_counts(
    intervals_path: 'Path',
    cram_path: 'CramPath',
    job_attrs: dict[str, str],
    output_base_path: str,
) -> 'BashJob':
    """
    Collect read counts from a CRAM file using GATK's CollectReadCounts tool, and the previously generated intervals

    Args:
        intervals_path (Path): path to the intervals file, generated in the previous step
        cram_path (CramPath): a SG's CRAM, Index, and Reference
        job_attrs (dict): collection of attributes to attach to the job
        output_base_path (str): root (minus file extensions) to write the output files to

    Returns:
        job (BashJob): The job object for the CollectReadCounts step
    """
    job = get_batch().new_bash_job('Collect gCNV read counts', job_attrs | {'tool': 'gatk CollectReadCounts'})
    job.image(image_path('gatk_gcnv'))

    reference = fasta_res_group(get_batch())

    # the specific GATK module here needs the full suffix ".counts.tsv.gz" to recognize the file format
    job.declare_resource_group(
        counts={
            'counts.tsv.gz': '{root}.counts.tsv.gz',
            'counts.tsv.gz.tbi': '{root}.counts.tsv.gz.tbi',
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
