from functools import cache
from typing import TYPE_CHECKING

from cpg_utils.cloud import read_secret
from cpg_utils.config import ConfigError, config_retrieve
from cpg_utils.hail_batch import get_batch
from google.api_core.exceptions import PermissionDenied
from loguru import logger

from cpg_gcnv.scripts import mt_to_es

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


@cache
def es_password() -> str:
    """
    Get Elasticsearch password. Moved into a separate method to simplify
    mocking in tests.
    """
    return read_secret(
        project_id=config_retrieve(['elasticsearch', 'password_project_id']),
        secret_name=config_retrieve(['elasticsearch', 'password_secret_id']),
        fail_gracefully=False,
    )


def submit_es_job_for_dataset(
    mt_path: str,
    index_name: str,
    done_flag: str,
    dataset: str,
) -> 'BashJob | None':
    """

    Args:
        mt_path ():
        index_name ():
        done_flag ():
        dataset ():

    Returns:

    """

    # try to generate a password here - we'll find out inside the script anyway, but
    # by that point we'd already have localised the MT, wasting time and money
    try:
        _es_password_string = es_password()
    except PermissionDenied:
        logger.warning(f'No permission to access ES password, skipping for {dataset}')
        return None
    except (ConfigError, KeyError):
        logger.warning(f'ES section not in config, skipping for {dataset}')
        return None

    # and just the name, used after localisation
    mt_name = mt_path.split('/')[-1]

    job = get_batch().new_bash_job(f'Generate {index_name} from {mt_path}')

    # set all job attributes in one bash
    job.cpu(4).memory('lowmem').storage('10Gi').image(config_retrieve(['workflow', 'driver_image']))

    # localise the MT
    job.command(f'gcloud --no-user-output-enabled storage cp -r {mt_path} $BATCH_TMPDIR')

    # run the export from the localised MT - this job writes no new data, just transforms and exports over network
    job.command(
        f"""
        python3 {mt_to_es.__file__} \
            --mt_path "${{BATCH_TMPDIR}}/{mt_name}" \
            --index {index_name!s} \
            --flag {done_flag!s}
        """,
    )

    return job
