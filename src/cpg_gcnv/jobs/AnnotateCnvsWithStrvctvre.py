from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def annotate_cnvs_with_strvctvre(
    input_vcf: str,
    output_vcf: str,
    job_attrs: dict[str, str],
) -> 'BashJob':
    """
    Annotate CNVs with STRVCTRE
    """

    job = hail_batch.get_batch().new_job('StrVCTVRE', job_attrs)

    job.image(config.config_retrieve(['images', 'strvctvre']))
    job.cpu(config.config_retrieve(['strvctvre_resources', 'cpu'], 2))
    job.memory(config.config_retrieve(['strvctvre_resources', 'memory'], '20Gi'))
    job.storage(config.config_retrieve(['strvctvre_resources', 'storage'], '20Gi'))

    strvctvre_phylop = config.config_retrieve(['annotate', 'strvctvre_phylop'])
    phylop_in_batch = hail_batch.get_batch().read_input(strvctvre_phylop)

    # read vcf and index into the batch
    input_vcf = hail_batch.get_batch().read_input_group(
        vcf=input_vcf,
        vcf_index=f'{input_vcf}.tbi',
    )['vcf']

    job.declare_resource_group(output_vcf={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    # run strvctvre
    job.command(f'python StrVCTVRE.py -i {input_vcf} -o temp.vcf -f vcf -p {phylop_in_batch}')
    job.command(f'bgzip temp.vcf -c > {job.output_vcf["vcf.bgz"]}')
    job.command(f'tabix {job.output_vcf["vcf.bgz"]}')

    hail_batch.get_batch().write_output(
        job.output_vcf,
        output_vcf.replace('.vcf.bgz', ''),
    )

    return job
