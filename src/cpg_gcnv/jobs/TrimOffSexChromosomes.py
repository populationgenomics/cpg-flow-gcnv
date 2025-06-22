from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def trim_sex_chromosomes(
    sgid_to_output: dict[str, Path],
    segment_vcfs: dict[str, Path],
    job_attrs: dict[str, str],
) -> 'list[BashJob]':
    """
    Read in a VCF, and remove all calls on X & Y.
    Relatively short-term solution to an ongoing gCNV issue, where gCNV's joint segmentation algorithm has been falsely
    detecting aneuploidies based on the inferred counts, and causing a halt to the workflow.
    An unresolved issue has been raised here https://github.com/broadinstitute/gatk/issues/8834
    """
    jobs = []

    # iterate over each of the SG/files we need to process
    for sgid, no_xy_vcf in sgid_to_output.items():
        if sgid == 'placeholder':
            continue

        sg_vcf = segment_vcfs[sgid]

        # create a new job for each SG
        job = hail_batch.get_batch().new_bash_job(
            f'Remove sex chromosomes from {sgid}',
            attributes=job_attrs | {'tool': 'bcftools'},
        )
        job.image(config.config_retrieve(['images', 'bcftools_120']))
        job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})
        localised_vcf = hail_batch.get_batch().read_input_group(
            **{
                'vcf.gz': sg_vcf,
                'vcf.gz.tbi': f'{sg_vcf}.tbi',
            },
        )['vcf.gz']
        autosomes = ' '.join([f'chr{i}' for i in range(1, 23)])
        job.command('set -euo pipefail')

        job.command(f'bcftools view -W=tbi -Oz -o {job.output["vcf.bgz"]} {localised_vcf} {autosomes}')
        hail_batch.get_batch().write_output(job.output, no_xy_vcf.removesuffix('.vcf.bgz'))
        jobs.append(job)

    return jobs
