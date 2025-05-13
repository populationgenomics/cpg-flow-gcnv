"""
Hail Query functions for seqr loader; CNV edition.
"""

import datetime
import gzip
from argparse import ArgumentParser
from os.path import join

from cpg_utils.hail_batch import genome_build, init_batch
from loguru import logger

import hail as hl

# I'm just going to go ahead and steal these constants from their seqr loader
GENE_SYMBOL = 'gene_symbol'
GENE_ID = 'gene_id'
MAJOR_CONSEQUENCE = 'major_consequence'

# Used to filter mt.info fields.
CONSEQ_PREDICTED_PREFIX = 'PREDICTED_'
NON_GENE_PREDICTIONS = {
    'PREDICTED_INTERGENIC',
    'PREDICTED_NONCODING_BREAKPOINT',
    'PREDICTED_NONCODING_SPAN',
}

GENCODE_FILE_HEADER = [
    'chrom',
    'source',
    'feature_type',
    'start',
    'end',
    'score',
    'strand',
    'phase',
    'info',
]


# yoinking some methods out of hail_scripts.computed_fields
# removes dependency on submodule completely
def get_expr_for_contig_number(locus: hl.LocusExpression) -> hl.Int32Expression:
    """Convert contig name to contig number"""
    return hl.bind(
        lambda contig: (
            hl.case().when(contig == 'X', 23).when(contig == 'Y', 24).when(contig[0] == 'M', 25).default(hl.int(contig))
        ),
        locus.contig.replace('^chr', ''),
    )


def get_expr_for_xpos(
    locus: hl.LocusExpression | hl.StructExpression,
) -> hl.Int64Expression:
    """Genomic position represented as a single number = contig_number * 10**9 + position.
    This represents chrom:pos more compactly and allows for easier sorting.
    """
    contig_number = get_expr_for_contig_number(locus)
    return hl.int64(contig_number) * 1_000_000_000 + locus.position


def parse_gtf_from_local(gtf_path: str) -> hl.dict:
    """
    Read over the localised GTF file and read into a dict

    n.b. due to a limit in Spark of 20MB per String length, the dictionary here is actually too large to be used
    in annotation expressions. To remedy this, the dictionary is returned as a list of fragments, and we can use each
    one in turn, then create a checkpoint between them.
    Args:
        gtf_path ():
        chunk_size (int): if specified, returns this dict as a list of dicts
    Returns:
        the gene lookup dictionary as a Hail DictExpression
    """
    gene_id_mapping = {}
    logger.info(f'Loading {gtf_path}')
    with gzip.open(gtf_path, 'rt') as gencode_file:
        # iterate over this file and do all the things
        for i, line in enumerate(gencode_file):
            line = line.rstrip('\r\n')
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) != len(GENCODE_FILE_HEADER):
                raise ValueError(f'Unexpected number of fields on line #{i}: {fields}')
            record = dict(zip(GENCODE_FILE_HEADER, fields, strict=True))
            if record['feature_type'] != 'gene':
                continue
            # parse info field
            info_fields_list = [x.strip().split() for x in record['info'].split(';') if x != '']
            info_fields = {k: v.strip('"') for k, v in info_fields_list}

            # skip an ENSG: ENSG mapping, redundant...
            if info_fields['gene_name'].startswith('ENSG'):
                continue
            gene_id_mapping[info_fields['gene_name']] = info_fields['gene_id'].split('.')[0]

    logger.info(f'Completed ingestion of gene-ID mapping, {len(gene_id_mapping)} entries')

    return [hl.literal(gene_id_mapping)]


def annotate_cohort_gcnv(vcf: str, mt_out: str, gencode: str, checkpoint: str):
    """
    Translate an annotated gCNV VCF into a Seqr-ready format
    Relevant gCNV specific schema
    https://github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/lib/model/gcnv_mt_schema.py
    Relevant gCNV loader script
    https://github.com/populationgenomics/seqr-loading-pipelines/blob/master/luigi_pipeline/seqr_gcnv_loading.py
    Args:
        vcf (str): Where is the VCF??
        mt_out (str): And where do you need output!?
        gencode (str): The path to a compressed GENCODE GTF file
        checkpoint (str): location we can write checkpoints to
    """

    init_batch()

    logger.info(f'Importing SV VCF {vcf}')
    mt = hl.import_vcf(
        vcf,
        array_elements_required=False,
        force_bgz=True,
        reference_genome=genome_build(),
        skip_invalid_loci=True,
    )

    # add attributes required for Seqr
    mt = mt.annotate_globals(
        sourceFilePath=vcf,
        genomeVersion=genome_build().replace('GRCh', ''),
        hail_version=hl.version(),
        datasetType='SV',
        sampleType='WES',
    )

    # reimplementation of
    # github.com/populationgenomics/seqr-loading-pipelines..luigi_pipeline/lib/model/gcnv_mt_schema.py
    mt = mt.annotate_rows(
        contig=mt.locus.contig.replace('^chr', ''),
        start=mt.locus.position,
        pos=mt.locus.position,
        # todo @MattWellie - review use of AC_Orig vs. AC (post-qc)
        sc=mt.info.AC_Orig[0],
        sf=mt.info.AF_Orig[0],
        sn=mt.info.AN_Orig,
        end=mt.info.END,
        sv_callset_Het=mt.info.N_HET,
        sv_callset_Hom=mt.info.N_HOMALT,
        gnomad_svs_ID=mt.info['gnomad_v2.1_sv_SVID'],
        gnomad_svs_AF=mt.info['gnomad_v2.1_sv_AF'],
        gnomad_svs_AC=hl.missing('float64'),
        gnomad_svs_AN=hl.missing('float64'),
        StrVCTVRE_score=hl.parse_float(mt.info.StrVCTVRE),
        svType=mt.alleles[1].replace('[<>]', ''),
        xpos=get_expr_for_xpos(mt.locus),
        xstart=get_expr_for_xpos(mt.locus),
        xstop=get_expr_for_xpos(hl.struct(contig=mt.locus.contig, position=mt.info.END)),
        num_exon=hl.agg.max(mt.NP),
    )

    # find all the column names which contain Gene symbols
    conseq_predicted_gene_cols = [
        gene_col
        for gene_col in mt.info
        if (gene_col.startswith(CONSEQ_PREDICTED_PREFIX) and gene_col not in NON_GENE_PREDICTIONS)
    ]

    # bank all those symbols before overwriting them - may not be required
    # might have to drop this for the Seqr ingest
    mt = mt.annotate_rows(
        geneSymbols=hl.set(
            hl.filter(
                hl.is_defined,
                [mt.info[gene_col] for gene_col in conseq_predicted_gene_cols],
            ).flatmap(
                lambda x: x,
            ),
        ),
    )

    mt = mt.checkpoint(join(checkpoint, 'pre-gene_annotation.mt'))

    # this next section is currently failing - the dictionary of genes is too large
    # to be used in an annotation expression. At least... I think it is
    # for i, chunks in enumerate(chunks(gene_id_mapping, 100)):

    # get the Gene-Symbol mapping dict
    gene_id_mapping = parse_gtf_from_local(gencode)

    # overwrite symbols with ENSG IDs in these columns
    # not sure why this is required, I think SV annotation came out
    # with ENSGs from the jump, but this is all symbols
    logger.info('Processing gene ID mapping chunk')
    for col_name in conseq_predicted_gene_cols:
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                **{col_name: hl.map(lambda gene: gene_id_mapping.get(gene, gene), mt.info[col_name])},
            ),
        )

    mt = mt.annotate_rows(
        # this expected mt.variant_name to be present, and it's not
        variantId=hl.format(f'%s_%s_{datetime.date.today():%m%d%Y}', mt.rsid, mt.svType),
        geneIds=hl.set(
            hl.filter(
                hl.is_defined,
                [mt.info[gene_col] for gene_col in conseq_predicted_gene_cols],
            ).flatmap(
                lambda x: x,
            ),
        ),
    )

    lof_genes = hl.set(mt.info.PREDICTED_LOF)
    major_consequence_genes = lof_genes | hl.set(mt.info.PREDICTED_COPY_GAIN)
    mt = mt.annotate_rows(
        sortedTranscriptConsequences=hl.map(
            lambda gene: hl.Struct(
                gene_id=gene,
                major_consequence=hl.or_missing(
                    major_consequence_genes.contains(gene),
                    hl.if_else(
                        lof_genes.contains(gene),
                        'LOF',
                        'COPY_GAIN',
                    ),
                ),
            ),
            mt.geneIds,
        ),
    )

    # transcriptConsequenceTerms
    default_consequences = [hl.format('gCNV_%s', mt.svType)]
    gene_major_consequences = hl.array(
        hl.set(
            mt.sortedTranscriptConsequences.filter(lambda x: hl.is_defined(x.major_consequence)).map(
                lambda x: x.major_consequence,
            ),
        ),
    )
    mt = mt.annotate_rows(
        transcriptConsequenceTerms=gene_major_consequences.extend(default_consequences),
        docId=mt.variantId[0:512],
    )

    # write this output
    mt.write(mt_out, overwrite=True)


def cli_main():
    """
    command line entrypoint
    """

    parser = ArgumentParser()
    parser.add_argument('--mt_out', help='Path to the MatrixTable, input or output', required=True)
    parser.add_argument('--checkpoint', help='Dir to write checkpoints to', required=True)
    parser.add_argument('--vcf', help='Path to input VCF file')
    parser.add_argument('--gencode', help='Path to input gencode GTF file')

    args = parser.parse_args()

    annotate_cohort_gcnv(vcf=args.vcf, mt_out=args.mt_out, gencode=args.gencode, checkpoint=args.checkpoint)


if __name__ == '__main__':
    cli_main()
