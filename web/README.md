# CCC Platform — web application

FastAPI back end (`backend/`), HTML pages with Alpine.js and plain CSS
(`frontend/`), and command-line tools (`tools/`). The scientific modules in
`../utils`, the pipeline scripts in `../scripts` and the knowledge bases in
`../data` are used as they are.

## Structure

```
web/
  run.py                        start the app (http://127.0.0.1:8000)
  backend/
    main.py                     routes: pages and /api/...
    config.py                   paths and settings
    datasets.py                 dataset registry and the active dataset
    dataset_admin.py            inspecting, preparing and analysing datasets
    jobs.py                     background jobs (preparation, CCI methods)
    settings.py                 analysis environment and Rscript detection
    index_store.py              fast SQLite index of the results
    score_definitions.py        what each tool's raw score means
    documentation.py            header and biology (UniProt, Reactome, SIGNOR)
    literature.py               PubMed search and evidence classification
    relevance.py                comparison of publications with the experiment
    expression.py               expression in the user's dataset (numpy only)
    project_context.py          experiment description (per dataset)
    assessments.py              expert assessments (../database/db.py)
    review_queue.py             review queue
    method_comparison.py        which tools report which interactions
    network.py                  communication network and matrix
    receptor_explorer.py        receptor recurrence (Fisher test, BH FDR)
    provenance.py               analysis records and draft Methods text
    mechanism.py                biological picture of an interaction (type, receptor, signalling)
    reactome_reactions.py       curated Reactome reactions of a ligand–receptor pair
  frontend/
    templates/                  pages; _doc_panel.html is shared
    static/css/app.css          design tokens and styles
    static/js/                  one file per page; common.js and mechanism.js are shared
  tools/
    prepare_dataset.py          check a .h5ad and store its expression matrix
    run_methods.py              run CellPhoneDB, NATMI, CellChat and the consensus
    rlauncher.py                start Rscript natively or inside WSL
    update_reactome_reactions.py  download Reactome's reaction mapping (once, internet)
  workspace/                    created at run time (index, settings, jobs)
```

Interactive API documentation: http://127.0.0.1:8000/api/docs

## Principles

* Each tool is used with its own criteria. Raw scores and ranks of different
  tools are never combined or compared; the platform shows which tools report
  an interaction.
* A computational prediction is a candidate, not evidence of binding or
  signalling. Knowledge about the specific ligand–receptor pair is kept apart
  from knowledge about the receptor in general, and from evidence in the
  user's own data.
* Complexes keep all their subunits throughout.

## Datasets

Each dataset is a folder `datasets/<id>/` with a `dataset.json`. Add one on the
Datasets page (or `python web/tools/prepare_dataset.py --help`), then run the
methods there (or `python web/tools/run_methods.py --id <id> --rscript <Rscript>`).
Each run writes `datasets/<id>/run/results/provenance.json`: software and
database versions, parameters, checksums of inputs and scripts, and timings.

## Settings (environment variables)

| Variable | Purpose |
|---|---|
| `NCBI_EMAIL` | Email sent to NCBI with PubMed requests (needed for new searches) |
| `CCC_ANALYSIS_PYTHON` | Python of the analysis environment (set by start.sh) |
| `CCC_RSCRIPT` | Rscript to use; `wsl:<path>` for an R inside WSL |
| `CCC_PUBMED_MAX_RESULTS` | Publications retrieved per pair (default 10) |
| `PORT` | Port (default 8000) |
