from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

if TYPE_CHECKING:
    from cpg_utils import Path


def counts_input_args(counts_paths: list['Path']) -> str:
    args = ''
    for f in counts_paths:
        counts = get_batch().read_input_group(
            **{
                'counts.tsv.gz': str(f),
                'counts.tsv.gz.tbi': str(f) + '.tbi',
            },
        )['counts.tsv.gz']
        args += f' --input {counts}'

    return args


def shard_items(name_only: bool):
    if shard_partition := config_retrieve(['workflow', 'interval_shards']):
        count = len(shard_partition)
        for i, interval in enumerate(shard_partition, start=1):
            name = 'part' + 'c'.join([s.removeprefix('chr') for s in interval])
            if name_only:
                yield name
            else:
                chroms = '|'.join(interval)
                yield name, i, count, f"awk '/^@/ || $1 ~ /^({chroms})$/'"
    else:
        raise NotImplementedError('gCNV sharding technique not specified')

