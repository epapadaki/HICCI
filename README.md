# CCC Platform

Cell–cell communication (CCC) analysis of single-cell RNA-seq data, with
biological documentation of every candidate interaction and expert review
(Human in the Loop).

The platform runs **CellPhoneDB**, **NATMI** and **CellChat**, each with its own
criteria, and shows which ligand–receptor interactions each tool reports. For
any interaction it documents the known biology of the ligand and receptor
(UniProt, Reactome, SIGNOR), the literature on the pair (PubMed), how close
that literature is to the user's experiment, and the expression of the partners
in the user's own data. Experts record their assessment of each interaction.

## Install (once)

Everything (web app, Python, R, CellPhoneDB, NATMI, CellChat) is installed into
one conda environment, `ccc-platform`.

* **Windows:** in PowerShell, from this folder:
  `wsl --cd "<this folder>" bash ./install.sh`
  (installs inside WSL/Ubuntu; WSL is included in Windows 10/11 — if missing,
  run `wsl --install` as administrator and restart).
* **Linux / macOS:** `./install.sh`

The first installation takes about 15–40 minutes (packages are downloaded and
CellChat is compiled). Conda (Miniforge) is installed automatically if absent.

## Start

* **Windows:** `wsl --cd "<this folder>" bash ./start.sh`
* **Linux / macOS:** `./start.sh`

Then open http://127.0.0.1:8000. Stop with Ctrl+C.

## Pages

| Page | What it is for |
|---|---|
| Datasets | Add a dataset (.h5ad), check it, run the CCI methods, see the analysis record |
| Interactions | All predictions, each tool with its own score and rank; biological documentation and expert assessment of any interaction |
| Network | Communication between cell populations (network and matrix views) |
| Receptors | Receptors that recur in a chosen set of interactions, with their ligands and populations |
| Review queue | Interactions grouped by which tools report them, and assessments that need a second look |
| Methods | Which interactions each tool reports, and their overlap |

## Example dataset

`datasets/endometrium_example/` describes the endometrium use case (the
CellPhoneDB example data). Put `endometrium_example_counts.h5ad` and
`endometrium_example_meta.tsv` in `data/cellphonedb/` to use it.

## Web interface only

To browse existing results without the analysis environment:
`pip install -r requirements.txt` and `python web/run.py`.

More detail on the code: `web/README.md`.
