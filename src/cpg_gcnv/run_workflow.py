#!/usr/bin/env python3

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow
from cpg_utils.config import config_retrieve

from cpg_gcnv.stages import MtToEsCnv, SplitAnnotatedCnvVcfByDataset


def cli_main():
    """
    CLI entrypoint - starts up the workflow
    """
    if config_retrieve(['workflow', 'sequencing_type']) != 'exome':
        raise ValueError('This is an exome-only workflow, for Genomes use GATK-SV instead')

    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    run_workflow(name='gcnv', stages=[MtToEsCnv, SplitAnnotatedCnvVcfByDataset], dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
