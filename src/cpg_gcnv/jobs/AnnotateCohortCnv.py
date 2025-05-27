from cpg_utils import config, hail_batch, Path

from cpg_gcnv.scripts import annotate_cohort


def submit_annotate_cohort_job(
    input_vcf: str,
    output_mt: str,
    checkpoint: Path,
    attributes: dict[str, str],
):
    """
    Submit a job to annotate a cohort with GATK SV
    """

    job = hail_batch.get_batch().new_job('AnnotateCohort gCNV multicohort', attributes)
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    gencode_gz = config.config_retrieve(['workflow', 'gencode_gtf_file'])
    gencode_gtf_local = hail_batch.get_batch().read_input(gencode_gz)

    job.command(
        f'python3 {annotate_cohort.__file__} '
        f'--mt_out {output_mt} '
        f'--checkpoint {checkpoint!s} '
        f'--vcf {input_vcf} '
        f'--gencode {gencode_gtf_local!s} ',
    )
    return job
