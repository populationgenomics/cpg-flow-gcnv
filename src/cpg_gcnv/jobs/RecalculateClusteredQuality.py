from typing import TYPE_CHECKING

from cpg_gcnv.utils import postprocess_calls

if TYPE_CHECKING:
    from cpg_utils import Path
    from hailtop.batch.job import Job


def recalculate_clustered_calls(
    ploidy_calls_path: 'Path',
    shard_paths: 'dict[str, Path]',
    sample_index: int,
    job_attrs: dict[str, str],
    output_prefix: str,
    clustered_vcf: str,
    intervals_vcf: str,
    qc_file: str,
) -> 'Job':
    """

    Args:
        ploidy_calls_path ():
        shard_paths ():
        sample_index ():
        job_attrs ():
        output_prefix ():
        clustered_vcf ():
        intervals_vcf ():
        qc_file ():

    Returns:

    """

    return postprocess_calls(
        ploidy_calls_path=ploidy_calls_path,
        shard_paths=shard_paths,
        sample_index=sample_index,
        job_attrs=job_attrs,
        output_prefix=output_prefix,
        clustered_vcf=clustered_vcf,
        intervals_vcf=intervals_vcf,
        qc_file=qc_file,
    )
