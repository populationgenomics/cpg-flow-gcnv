import functools
from typing import TYPE_CHECKING

import loguru
from cpg_utils import cloud, config, hail_batch
from google.api_core.exceptions import PermissionDenied

from cpg_gcnv.scripts import mt_to_es

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


@functools.cache
def es_password() -> str:
    """
    Get Elasticsearch password. Moved into a separate method to simplify
    mocking in tests.
    """
    return cloud.read_secret(
        project_id=config.config_retrieve(['elasticsearch', 'password_project_id']),
        secret_name=config.config_retrieve(['elasticsearch', 'password_secret_id']),
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
        loguru.logger.warning(f'No permission to access ES password, skipping for {dataset}')
        return None
    except (config.ConfigError, KeyError):
        loguru.logger.warning(f'ES section not in config, skipping for {dataset}')
        return None

    # and just the name, used after localisation
    mt_name = mt_path.split('/')[-1]

    job = hail_batch.get_batch().new_bash_job(f'Generate {index_name} from {mt_path}')

    # set all job attributes in one bash
    job.cpu(4).memory('lowmem').storage('10Gi').image(config.config_retrieve(['workflow', 'driver_image']))

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
