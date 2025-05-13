from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve, image_path
from cpg_utils.hail_batch import get_batch

from cpg_gcnv.scripts import update_vcf_attributes

if TYPE_CHECKING:
    from cpg_utils import Path
    from hailtop.batch.job import BashJob


def fast_merge_calls(sg_vcfs: list[str], job_attrs: dict[str, str], output_path: 'Path') -> list['BashJob']:
    """
    This job will run a fast simple merge on per-SGID call files
    It then throws in a python script to add in two additional header lines
    and edit the SVLEN and SVTYPE attributes into each row

    Escape here for single-VCF merges - the merge command is invalid, so we pass through the first step

    Args:
        sg_vcfs (list[str]): paths to all individual VCFs
        job_attrs (dict): any params to attach to the job
        output_path (Path): path to the final merged VCF
    """

    merge_job = get_batch().new_job('Merge gCNV calls', job_attrs | {'tool': 'bcftools'})
    merge_job.image(image_path('bcftools_120'))

    # this should be made reactive, in case we scale past 10GB
    merge_job.storage('10Gi')

    batch_vcfs = []
    for each_vcf in sg_vcfs:
        batch_vcfs.append(
            get_batch().read_input_group(**{'vcf.gz': each_vcf, 'vcf.gz.tbi': f'{each_vcf}.tbi'})['vcf.gz'],
        )

    if len(batch_vcfs) == 0:
        raise ValueError('No VCFs to merge')

    if len(batch_vcfs) == 1:
        merge_job.tmp_vcf = batch_vcfs[0]

    else:
        # option breakdown:
        # -Oz: bgzip output
        # -o: output file
        # --threads: number of threads to use
        # -m: merge strategy
        # -0: compression level
        merge_job.command(f'bcftools merge {" ".join(batch_vcfs)} -Oz -o {merge_job.tmp_vcf} --threads 4 -m all -0')

    # now normalise the result, splitting multiallelics
    merge_job.command(f'bcftools norm -m -any -Oz -o {merge_job.tmp_vcf_split} {merge_job.tmp_vcf}')

    # second job, do the VCF content updates
    update_job = get_batch().new_bash_job('Update VCF content')
    update_job.image(config_retrieve(['workflow', 'driver_image']))
    update_job.storage('10Gi')
    update_job.command(f'{update_vcf_attributes.__file__} {merge_job.tmp_vcf_split} {update_job.output}')

    # a third job just to tidy up
    third_job = get_batch().new_job('bgzip and tabix')
    third_job.image(image_path('bcftools_120'))
    third_job.storage('10Gi')
    third_job.declare_resource_group(output={'vcf.bgz': '{root}.vcf.bgz', 'vcf.bgz.tbi': '{root}.vcf.bgz.tbi'})

    third_job.command(f'bcftools view -W=tbi -Oz -o {third_job.output["vcf.bgz"]} {update_job.output}')

    # get the output root to write to
    get_batch().write_output(third_job.output, str(output_path).removesuffix('.vcf.bgz'))
    return [merge_job, update_job, third_job]
