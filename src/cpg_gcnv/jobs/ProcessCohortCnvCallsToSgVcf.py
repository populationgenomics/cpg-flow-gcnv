from typing import TYPE_CHECKING

from cpg_utils import Path

from cpg_gcnv.utils import postprocess_calls

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def postprocess_unclustered_calls(
    ploidy_calls_path: Path,
    shard_paths: dict[str, Path],
    sample_index: int,
    job_attrs: dict[str, str],
    output_prefix: str,
) -> 'BashJob':
    """

    Args:
        ploidy_calls_path ():
        shard_paths ():
        sample_index ():
        job_attrs ():
        output_prefix ():
    """

    return postprocess_calls(
        ploidy_calls_path=ploidy_calls_path,
        shard_paths=shard_paths,
        sample_index=sample_index,
        job_attrs=job_attrs,
        output_prefix=output_prefix,
    )
