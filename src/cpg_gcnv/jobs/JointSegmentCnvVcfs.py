from typing import TYPE_CHECKING

from cpg_flow import resources
from cpg_utils import Path, config, hail_batch

from cpg_gcnv.utils import chunks

if TYPE_CHECKING:
    from hailtop.batch import Resource, ResourceFile, ResourceGroup
    from hailtop.batch.job import BashJob


def joint_segment_vcfs(
    segment_vcfs: list['ResourceFile'],
    pedigree: 'ResourceFile',
    reference: 'ResourceGroup',
    intervals: 'ResourceFile',
    title: str,
    job_attrs: dict,
) -> tuple['BashJob', 'Resource']:
    """
    This job will run the joint segmentation step of the gCNV workflow
    Takes individual Segment VCFs and merges them into a single VCF
    Depending on the config setting workflow.num_samples_per_scatter_block
    this may be conducted in hierarchical 2-step, with intermediate merges
    being conducted, then a merge of those intermediates

    Returns:
        the job that does the work, and the resulting resource group of VCF & index
    """
    job = hail_batch.get_batch().new_bash_job(f'Joint Segmentation {title}', job_attrs | {'tool': 'gatk'})
    job.declare_resource_group(output={'vcf.gz': '{root}.vcf.gz', 'vcf.gz.tbi': '{root}.vcf.gz.tbi'})
    job.image(config.image_path('gatk_gcnv'))

    # set highmem resources for this job
    job_res = resources.HIGHMEM.request_resources(ncpu=2, storage_gb=10)
    job_res.set_to_job(job)

    vcf_string = ''
    for each_vcf in segment_vcfs:
        vcf_string += f' -V {each_vcf}'

    # this already creates a tabix index
    job.command(
        f"""
    set -e
    gatk --java-options "{job_res.java_mem_options()}" JointGermlineCNVSegmentation \\
        -R {reference.base} \\
        -O {job.output['vcf.gz']} \\
        {vcf_string} \\
        --model-call-intervals {intervals} \\
        -ped {pedigree}
    """,
    )
    return job, job.output


def run_joint_segmentation(
    segment_vcfs: list[str],
    pedigree: Path,
    intervals: Path,
    tmp_prefix: Path,
    output_path: Path,
    job_attrs: dict[str, str] | None = None,
) -> 'list[BashJob]':
    """
    This job will run the joint segmentation step of the gCNV workflow
    Takes individual Segment VCFs and merges them into a single VCF
    Depending on the config setting workflow.num_samples_per_scatter_block
    this may be conducted in hierarchical 2-step, with intermediate merges
    being conducted, then a merge of those intermediates

    Args:
        segment_vcfs ():
        pedigree ():
        intervals ():
        tmp_prefix ():
        output_path ():
        job_attrs ():

    Returns:

    """
    jobs = []

    pedigree_in_batch = hail_batch.get_batch().read_input(pedigree)
    intervals_in_batch = hail_batch.get_batch().read_input(intervals)

    # find the number of samples to shove into each scatter block
    sams_per_block = config.config_retrieve(['workflow', 'num_samples_per_scatter_block'])

    reference = hail_batch.fasta_res_group(hail_batch.get_batch())

    chunked_vcfs = []

    # if we have more samples to process than the block size, condense
    # this is done by calling an intermediate round of VCF segmenting
    if len(segment_vcfs) > sams_per_block:
        for subchunk_index, chunk_vcfs in enumerate(chunks(segment_vcfs, sams_per_block)):
            # create a new job for each chunk
            # read these files into this batch
            local_vcfs = [
                hail_batch.get_batch().read_input_group(
                    vcf=vcf,
                    index=f'{vcf}.tbi',
                )['vcf']
                for vcf in chunk_vcfs
            ]
            job, vcf_group = joint_segment_vcfs(
                local_vcfs,
                pedigree=pedigree_in_batch,
                reference=reference,
                intervals=intervals_in_batch,
                job_attrs=job_attrs or {} | {'title': f'sub-chunk_{subchunk_index}'},
                title=f'sub-chunk_{subchunk_index}',
            )
            chunked_vcfs.append(vcf_group['vcf.gz'])
            hail_batch.get_batch().write_output(vcf_group, tmp_prefix / f'subchunk_{subchunk_index}')
            jobs.append(job)

    # else, all vcf files into the batch
    else:
        chunked_vcfs = [
            hail_batch.get_batch()
            .read_input_group(
                vcf=vcf,
                index=f'{vcf}.tbi',
            )
            .vcf
            for vcf in segment_vcfs
        ]

    # second round of condensing output - produces one single file
    job, vcf_group = joint_segment_vcfs(
        chunked_vcfs,
        pedigree=pedigree_in_batch,
        reference=reference,
        intervals=intervals_in_batch,
        job_attrs=job_attrs or {} | {'title': 'all-chunks'},
        title='all-chunks',
    )
    jobs.append(job)

    # write the final output file group (VCF & index)
    hail_batch.get_batch().write_output(vcf_group, f'{str(output_path).removesuffix(".vcf.gz")}')
    return jobs
