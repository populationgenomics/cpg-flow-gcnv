"""
Hail Query functions for seqr loader; CNV edition.
"""

from argparse import ArgumentParser

from cpg_utils import hail_batch
import loguru

import hail as hl


def annotate_dataset_gcnv(mt_in: str, mt_out: str):
    """
    process data specific to samples in this dataset
    do this after sub-setting to specific samples
    Args:
        mt_in (str): path to the annotated MatrixTable
        mt_out (str): and where do you want it to end up?
    """

    hail_batch.init_batch()

    loguru.logger.info('Annotating genotypes')

    mt = hl.read_matrix_table(mt_in)

    # adding in the GT here, that may cause problems later?
    mt = mt.annotate_rows(
        genotypes=hl.agg.collect(
            hl.struct(
                sample_id=mt.s,
                gq=mt.QA,
                cn=mt.CN,
                end=mt.end,
                num_exon=mt.NP,
                start=mt.start,
                geneIds=mt.geneIds,
                gt=mt.GT,
            ),
        ),
    )

    def _genotype_filter_samples(fn) -> hl.set[str]:
        # Filter on the genotypes.
        return hl.set(mt.genotypes.filter(fn).map(lambda g: g.sample_id))

    def _genotype_filter_samples_cn2() -> hl.set[str]:
        # Filter on the genotypes.
        return hl.set(mt.genotypes.filter(lambda g: ((g.gt.is_haploid()) & (g.cn == 2))).map(lambda g: g.sample_id))  # noqa: PLR2004

    # top level - decorator
    def _capture_i_decorator(func):  # noqa: ANN202
        # call the returned_function(i) which locks in the value of i
        def _inner_filter(i):  # noqa: ANN202
            # the _genotype_filter_samples will call this _func with g
            def _func(g):  # noqa: ANN202
                return func(i, g)

            return _func

        return _inner_filter

    @_capture_i_decorator
    def _filter_sample_cn(i, g) -> bool:
        return g.cn == i

    @_capture_i_decorator
    def _filter_samples_gq(i, g) -> bool:
        return (g.gq >= i) & (g.gq < i + 10)

    @_capture_i_decorator
    def _filter_num_alt(i, g) -> bool:
        return i == g.gt.n_alt_alleles()

    # github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/gcnv_mt_schema.py
    mt = mt.annotate_rows(
        # dubious about this annotation - expected field is qs, I'm using gq, derived from CNQ
        samples_qs=hl.struct(
            **{f'{i}_to_{i + 10}': _genotype_filter_samples(_filter_samples_gq(i)) for i in range(0, 1000, 10)},
            gt_1000=_genotype_filter_samples(lambda g: g.gq >= 1000),  # noqa: PLR2004
        ),
        # ok, here's what we're
        samples_cn=hl.struct(
            **{str(i): _genotype_filter_samples(_filter_sample_cn(i)) for i in [0, 1, 3]},
            gte_4=_genotype_filter_samples(lambda g: g.cn >= 4),  # noqa: PLR2004
            # and a special case for male sex-chrom CN
            **{'2': _genotype_filter_samples_cn2()},
        ),
        # re-embedding the samples_num_alt, derived from hl.call().n_alt_alleles()
        samples_num_alt=hl.struct(**{str(i): _genotype_filter_samples(_filter_num_alt(i)) for i in range(1, 3, 1)}),
    )

    # removing GT again, out of an abundance of caution
    # adding in the GT here, that may cause problems later?
    mt = mt.annotate_rows(
        genotypes=hl.agg.collect(
            hl.struct(
                sample_id=mt.s,
                gq=mt.QA,
                cn=mt.CN,
                end=mt.end,
                num_exon=mt.NP,
                start=mt.start,
                geneIds=mt.geneIds,
            ),
        ),
    )
    loguru.logger.info('Genotype fields annotated')
    mt.describe()
    mt.write(mt_out, overwrite=True)
    loguru.logger.info(f'Written gCNV MT to {mt_out}')


def cli_main():
    """
    command line entrypoint
    """

    parser = ArgumentParser()
    parser.add_argument('--mt_in', help='Path to input MT', required=True)
    parser.add_argument('--mt_out', help='Path to the MatrixTable, input or output', required=True)
    args = parser.parse_args()

    annotate_dataset_gcnv(mt_in=args.mt_in, mt_out=args.mt_out)


if __name__ == '__main__':
    cli_main()
