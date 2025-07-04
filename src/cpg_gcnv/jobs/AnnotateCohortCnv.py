from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def submit_annotate_cohort_job(
    input_vcf: str,
    output_mt: str,
    checkpoint: Path,
    attributes: dict[str, str],
) -> 'BashJob':
    """Submit a job to annotate a cohort with GATK SV."""
    batch_instance = hail_batch.get_batch()

    job = batch_instance.new_job('AnnotateCohort gCNV multicohort', attributes)
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    gencode_gz = config.config_retrieve(['annotate', 'gencode_gtf_file'])
    gencode_gtf_local = batch_instance.read_input(gencode_gz)

    job.command(
        f"""
        python3 -m cpg_gcnv.scripts.annotate_cohort \\
        --mt_out {output_mt} \\
        --checkpoint {checkpoint!s} \\
        --vcf {input_vcf} \\
        --gencode {gencode_gtf_local!s}
        """,
    )
    return job
