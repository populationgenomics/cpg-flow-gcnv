#!/usr/bin/env python3

"""
This is the main entry point for the workflow.
This is a re-implementation of the canonical main.py file in production-pipelines.
The purpose of this script is to import all the Stages in the workflow (or at least the terminal workflow nodes)
and begin the CPG-Flow Stage discovery and graph construction process

This is re-implemented as a simpler form, only knowing how to build a single workflow, instead of choosing at runtime
"""

from argparse import ArgumentParser

from cpg_flow.workflow import run_workflow
from cpg_utils.config import config_retrieve

from cpg_gcnv.stages import CallGermlineCnvsWithGatk


def cli_main():
    """
    CLI entrypoint - starts up the workflow
    """
    assert config_retrieve(['workflow', 'sequencing_type']) == 'exome', 'For exomes, use GATK-SV instead'

    parser = ArgumentParser()
    parser.add_argument('--dry_run', action='store_true', help='Dry run')
    args = parser.parse_args()

    run_workflow(stages=[CallGermlineCnvsWithGatk], dry_run=args.dry_run)


if __name__ == '__main__':
    cli_main()
