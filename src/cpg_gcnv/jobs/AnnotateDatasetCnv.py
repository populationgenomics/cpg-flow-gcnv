from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch

from cpg_gcnv.scripts import annotate_dataset

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def submit_annotate_dataset_job(
    input_mt: str,
    output_mt: str,
    attributes: dict[str, str],
) -> 'BashJob':
    """
    Submit a job to annotate a cohort with GATK SV
    """

    job = hail_batch.get_batch().new_bash_job('AnnotateDataset gCNV multicohort', attributes)
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    job.command(f'python3 {annotate_dataset.__file__} --mt_in {input_mt} --mt_out {output_mt}')

    return job
