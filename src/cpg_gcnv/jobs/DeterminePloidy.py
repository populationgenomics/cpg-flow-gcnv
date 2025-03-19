from typing import TYPE_CHECKING

from cpg_utils.config import image_path
from cpg_utils.hail_batch import fasta_res_group, get_batch
from cpg_flow.resources import HIGHMEM

if TYPE_CHECKING:
    from cpg_utils import Path
    from cpg_flow.filetypes import CramPath
    from hailtop.batch.job import Job


def _counts_input_args(counts_paths: list[Path]) -> str:
    args = ''
    for f in counts_paths:
        counts = get_batch().read_input_group(
            **{
                'counts.tsv.gz': str(f),
                'counts.tsv.gz.tbi': str(f) + '.tbi',
            },
        )["counts.tsv.gz"]
        args += f' --input {counts}'

    return args

def filter_and_determine_ploidy(
    ploidy_priors_path: str,
    preprocessed_intervals_path: 'Path',
    annotated_intervals_path: 'Path',
    counts_paths: list['Path'],
    job_attrs: dict[str, str],
    output_paths: dict[str, 'Path'],
) -> list['Job']:
    j = get_batch().new_job(
        'Filter intervals and determine ploidy',
        job_attrs
        | {
            'tool': 'gatk FilterIntervals/DetermineGermlineContigPloidy',
        },
    )
    j.image(image_path('gatk_gcnv'))

    # set highmem resources for this job
    job_res = HIGHMEM.request_resources(ncpu=2, storage_gb=10)
    job_res.set_to_job(j)

    counts_input_args = _counts_input_args(counts_paths)
    cmd = ''

    if can_reuse(output_paths['filtered']):
        # Remove 'filtered' entry from output_paths so we don't write_output() it later
        filtered: ResourceFile = get_batch().read_input(str(output_paths.pop('filtered')))
    else:
        preprocessed_intervals = get_batch().read_input(str(preprocessed_intervals_path))
        annotated_intervals = get_batch().read_input(str(annotated_intervals_path))

        cmd += f"""
        gatk --java-options "{job_res.java_mem_options()}" FilterIntervals \\
          --interval-merging-rule OVERLAPPING_ONLY \\
          --intervals {preprocessed_intervals} --annotated-intervals {annotated_intervals} \\
          {counts_input_args} \\
          --output {j.filtered}
        """

        assert isinstance(j.filtered, JobResourceFile)
        j.filtered.add_extension('.interval_list')
        filtered = j.filtered

    # (Other arguments may be cloud URLs, but this *must* be a local file)
    ploidy_priors = get_batch().read_input(ploidy_priors_path)

    cmd += f"""
    gatk --java-options "{job_res.java_mem_options()}" DetermineGermlineContigPloidy \\
      --interval-merging-rule OVERLAPPING_ONLY \\
      --intervals {filtered} --contig-ploidy-priors {ploidy_priors} \\
      {counts_input_args} \\
      --output $BATCH_TMPDIR --output-prefix ploidy

    tar -czf {j.calls} -C $BATCH_TMPDIR ploidy-calls
    tar -czf {j.model} -C $BATCH_TMPDIR ploidy-model
    """

    j.command(command(cmd))
    for key, path in output_paths.items():
        get_batch().write_output(j[key], str(path))
    return [j]
