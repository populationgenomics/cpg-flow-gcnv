from typing import TYPE_CHECKING

from cpg_flow import targets, workflow
from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def split_mc_vcf_by_dataset(
    dataset: targets.Dataset,
    input_vcf: str,
    output: Path,
    job_attrs: dict[str, str],
) -> 'BashJob':
    """Use BCFtools to split a multicohort VCF by dataset."""

    # write a temp file for this dataset containing all relevant SGs
    sgids_list_path = dataset.tmp_prefix() / workflow.get_workflow().output_version / 'sgid-list.txt'
    if not config.config_retrieve(['workflow', 'dry_run'], False):
        with sgids_list_path.open('w') as f:
            for sgid in dataset.get_sequencing_group_ids():
                f.write(f'{sgid}\n')

    job = hail_batch.get_batch().new_job(
        name=f'SplitAnnotatedVcfByDataset: {dataset}',
        attributes=job_attrs | {'tool': 'bcftools'},
    )
    job.image(config.config_retrieve(['images', 'bcftools_120'])).cpu(1).memory('highmem').storage('10Gi')

    job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    local_sgid_file = hail_batch.get_batch().read_input(sgids_list_path)
    job.command(
        f"""
        bcftools view \\
            {input_vcf} \\
            --force-samples \\
            -S {local_sgid_file} \\
            -Oz \\
            -o {job.output['vcf.bgz']} \\
            -W=tbi
        """,
    )

    hail_batch.get_batch().write_output(job.output, str(output).removesuffix('.vcf.bgz'))

    return job
