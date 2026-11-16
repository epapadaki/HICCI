# HICCI — Human-In-the-loop Curation of Cellular Interactions

HICCI is a web platform for the expert curation of cell–cell interactions
predicted from single-cell RNA-seq data. It runs an established inference
method chosen by the user, ranks the candidate interactions within each method
on a percentile scale and pre-selects those above a user-defined cutoff. The
expert then reviews each candidate against separate lines of evidence, keeps,
removes or rescues it, and every decision is recorded. Network, pathway and
enrichment analyses are recomputed on the curated set after each decision.

## Main features

- **Inference methods:** CellPhoneDB v5, NATMI, CellChat v2, SingleCellSignalR (v1) and
  ICELLNET, run with their own databases and reporting criteria. One method is
  enough; several can be run on the same dataset.
- **Complex-aware interactions:** multi-subunit ligands and receptors are kept
  as complexes, and an interaction is defined by sender, receiver, ligand
  entity and receptor entity.
- **Within-method percentile ranking** and a default set defined by a cutoff
  chosen by the user.
- **Five separate lines of evidence for every candidate**, never merged into a
  single score:
  - prior knowledge (UniProt, Reactome, SIGNOR, CellChatDB);
  - literature support (PubMed);
  - relevance of the literature to the study;
  - expression in the user's data;
  - signalling response in the receiver (PROGENy, CollecTRI).
- **Decision record:** every keep, remove or rescue decision is stored and can
  be exported.
- **Downstream analyses on the curated set:**
  - curated network;
  - communication profiles;
  - Reactome processes;
  - molecular hubs;
  - ligand–receptor graphs;
  - functional enrichment.
- **Optional AI-assisted components (Google Gemini):**
  - method advisor;
  - literature classification;
  - Discussion draft.

  They support the expert but take no part in the curation decisions. Without
  an API key, the platform uses its rule-based behaviour.

## Installation

HICCI runs on Linux, macOS and Windows (through WSL for the R-based methods).
All components are installed in one conda environment.

1. Install [Miniforge](https://github.com/conda-forge/miniforge) (conda).
2. Clone the repository and create the environment:

   ```bash
   git clone https://github.com/epapadaki/HICCI.git
   cd HICCI
   conda env create -f install/environment.yml
   conda activate ccc-platform
   ```

3. Install the R packages that conda cannot provide (CellChat, NMF, presto, ICELLNET):

   ```bash
   Rscript install/install_r_packages.R
   Rscript install/install_icellnet.R
   ```

4. Install NATMI (only needed if NATMI is used):

   ```bash
   git clone https://github.com/forrest-lab/NATMI.git external/NATMI
   git -C external/NATMI checkout c7ad5a9
   ```

   NATMI is used unmodified. Its last version (2021) calls a pandas function
   that pandas 2 removed; HICCI handles this automatically when it runs NATMI
   (`web/tools/natmi_compat`), so no manual change is needed.

The CellPhoneDB database (v5.0.0) is downloaded automatically on the first
CellPhoneDB run. The first installation takes about 15–40 minutes, mainly
because CellChat is compiled.

### Optional settings

Copy `.env.example` to `.env` and fill in the values you need. The `.env` file
is ignored by git and is never published.

| Variable | Purpose |
|---|---|
| `GEMINI_API_KEY` | Enables the AI-assisted components |
| `GEMINI_MODEL` | Gemini model (default: `gemini-2.5-flash`) |
| `NCBI_EMAIL`, `NCBI_API_KEY` | Recommended by NCBI for PubMed queries |

## Start

```bash
conda activate ccc-platform
python web/run.py
```

Then open <http://127.0.0.1:8000> in a browser. Stop with Ctrl+C.

On Windows, `start.bat` or `start-ccc.ps1` starts the web interface from a
local Python environment (`.venv`).

## Workflow

| Page | Purpose |
|---|---|
| Datasets | Upload an annotated `.h5ad` file, select the cell-type column, check the data and describe the experiment |
| Tools | Choose and run the inference method(s), with the method advisor |
| Interactions | Ranked candidates, the evidence for each one, and the keep / remove / rescue decisions |
| Curated results | Network, communication profiles, processes, hubs and ligand–receptor views of the curated set |
| Functional enrichment | Over-representation analysis of the curated ligand and receptor genes |
| Discussion | AI-generated draft Discussion of the curated results (optional) |

## Example dataset (melanoma)

The case study of the paper uses the melanoma tumour microenvironment dataset
of Tirosh et al. (2016, GEO GSE72056; 4,139 cells, 7 populations).

1. Download the prepared file `melanoma_tme.h5ad` (295 MB) from Zenodo:
   <https://doi.org/10.5281/zenodo.23976292>
   (MD5: `4be06f07c5d5098a74007f4f2f9affa4`)
2. Place it at `data/uploads/melanoma_tme.h5ad`.
3. Start the platform, select the dataset on the Datasets page and run CellChat
   from the Tools page (about 40 minutes).

## Repository structure

```
web/        web application (FastAPI back end, HTML/Alpine.js front end, tools)
utils/      scientific modules (knowledge bases, literature, Gemini client)
scripts/    input preparation, method runners and output standardisation
database/   storage of expert decisions (SQLite)
data/       local snapshots of the knowledge bases and gene-set libraries
datasets/   dataset descriptions (results are written here at run time)
install/    conda environment and R installation scripts
```

## Citation

If you use HICCI, please cite:

> Papadaki E., Vrahatis A.G. HICCI: a human-in-the-loop platform for the expert
> curation of predicted cell–cell interactions. (submitted)

## License

MIT License — see [LICENSE](LICENSE).
