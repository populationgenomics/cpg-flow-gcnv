# gCNV, Using CPG-Flow

This is the gCNV workflow, migrated from [Production-Pipelines](https://github.com/populationgenomics/production-pipelines/blob/main/cpg_workflows/stages/gcnv.py) to the CPG-Flow framework.

Current Version: 0.1.3

## Purpose

When migrating workflows from production-pipelines, this template respository structure can be used to start with a
sensible directory structure, and some suggested conventions for naming and placement of files.

```txt
src
├── cpg_gcnv
│   ├── __init__.py
│   ├── config_template.toml
│   ├── jobs
│   │   ├── AnnotateCnvsWithStrvctvre.py
│   │   ├── AnnotateCnvsWithSvAnnotate.py
│   │   ├── AnnotateCohortCnv.py
│   │   ├── ...
│   ├── run_workflow.py
│   ├── scripts
│   │   ├── __init__.py
│   │   ├── annotate_cohort.py
│   │   ├── annotate_dataset.py
│   │   ├── mt_to_es.py
│   │   ├── ...
│   ├── stages.py
│   └── utils.py
```

`stages.py` contains Stages in the workflow, with the actual logic imported from files in `jobs`.

`jobs/` is a folder containing the logic for each of the Stages in the workflow. Each file contains a single Stage, and has a name mirroring the Stage it contains. This implements all the logic for a stage, including assembling the input dictionary to be passed to Cromwell and scheduling any Batch Jobs.

`config_template.toml` is a base config, indicating settings which are mandatory for the pipeline to run. Actual workflow configs should inherit from this, and can add additional settings as required (e.g. dataset, cohort, updated image paths)

Example Analysis-Runner invocation:

```bash
analysis-runner \
    --skip-repo-checkout \
    --image australia-southeast1-docker.pkg.dev/cpg-common/images/cpg-flow-gcnv:0.1.3 \
    --dataset DATASET \
    --description 'gCNV, CPG-flow' \
    -o gCNV_cpg-flow \
    --access-level full \
    --config CONFIG \
    run_workflow
```
