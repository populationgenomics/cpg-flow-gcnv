from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

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

    job = get_batch().new_bash_job('AnnotateDataset gCNV multicohort', attributes)
    job.image(config_retrieve(['workflow', 'driver_image']))

    job.command(
        f'{annotate_dataset.__file__} '
        f'--mt_in {input_mt} '
        f'--mt_out {output_mt} ',
    )

    return job
