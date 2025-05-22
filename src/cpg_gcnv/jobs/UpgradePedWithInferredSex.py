from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch

from cpg_gcnv.scripts import upgrade_ped_with_inferred_sex

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def upgrade_ped_file(
    ped_file: str,
    new_output: str,
    aneuploidies: str,
    ploidy_tar: str,
) -> 'BashJob':
    """
    Update the default Pedigree with the inferred ploidy information
    update a ped file, and

    Args:
        ped_file (str): path to a PED in GCP
        new_output ():
        aneuploidies (str): where to write identified aneuploidies
        ploidy_tar ():
    """

    local_ped = hail_batch.get_batch().read_input(ped_file)
    j = hail_batch.get_batch().new_bash_job('Upgrade PED file with inferred Ploidy')
    j.image(config.config_retrieve(['workflow', 'driver_image']))

    # path to the python script
    j.command(f'tar -xf {ploidy_tar} -C .')  # creates the folder ploidy-calls
    j.command(f'python3 {upgrade_ped_with_inferred_sex.__file__} {local_ped} {j.output} {j.aneuploidies} ploidy-calls')
    hail_batch.get_batch().write_output(j.output, new_output)
    hail_batch.get_batch().write_output(j.aneuploidies, aneuploidies)
    return j
