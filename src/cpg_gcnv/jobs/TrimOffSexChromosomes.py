from typing import TYPE_CHECKING

from cpg_utils.config import image_path
from cpg_utils.hail_batch import get_batch

if TYPE_CHECKING:
    from cpg_utils import Path
    from hailtop.batch.job import Job


def trim_sex_chromosomes(
    sgid_to_output: 'dict[str, Path]',
    segment_vcfs: 'dict[str, Path]',
    job_attrs: dict[str, str],
) -> 'list[Job]':
    """

    Args:
        sgid_to_output ():
        segment_vcfs ():
        job_attrs ():

    Returns:

    """
    jobs = []

    # iterate over each of the SG/files we need to process
    for sgid, no_xy_vcf in sgid_to_output.items():
        if sgid == 'placeholder':
            continue

        sg_vcf = segment_vcfs[sgid]

        # create a new job for each SG
        job = get_batch().new_bash_job(f'Remove sex chromosomes from {sgid}', job_attrs | {'tool': 'bcftools'})
        job.image(image_path('bcftools_120'))
        job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
        localised_vcf = get_batch().read_input_group(**{'vcf.gz': sg_vcf, 'vcf.gz.tbi': f'{sg_vcf}.tbi'})['vcf.gz']
        autosomes = ' '.join([f'chr{i}' for i in range(1, 23)])
        job.command('set -euo pipefail')

        # TODO when we adopt bcftools 1.20+ as standard we can drop the separate tabix step
        job.command(f'bcftools view -W=tbi -Oz -o {job.output["vcf.bgz"]} {localised_vcf} {autosomes}')
        get_batch().write_output(job.output, no_xy_vcf.removesuffix('.vcf.bgz'))
        jobs.append(job)

    return jobs
