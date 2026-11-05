from pathlib import Path

from cellphonedb.src.core.methods import cpdb_analysis_method
from cellphonedb.utils import db_utils


# --------------------------------------------------
# Project paths
# --------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[1]

DATA_DIR = PROJECT_ROOT / "data" / "cellphonedb"
RESULTS_DIR = PROJECT_ROOT / "results" / "cellphonedb"
DB_DIR = DATA_DIR / "database" / "v5.0.0"

COUNTS_FILE = DATA_DIR / "endometrium_example_counts.h5ad"
META_FILE = DATA_DIR / "endometrium_example_meta.tsv"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
DB_DIR.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# Check input files
# --------------------------------------------------
if not COUNTS_FILE.exists():
    raise FileNotFoundError(
        f"Counts file not found: {COUNTS_FILE}"
    )

if not META_FILE.exists():
    raise FileNotFoundError(
        f"Metadata file not found: {META_FILE}"
    )


# --------------------------------------------------
# Download CellPhoneDB database if needed
# --------------------------------------------------
CPDB_VERSION = "v5.0.0"
CPDB_FILE = DB_DIR / "cellphonedb.zip"

if not CPDB_FILE.exists():
    print("Downloading CellPhoneDB database...")
    db_utils.download_database(
        str(DB_DIR),
        CPDB_VERSION
    )

if not CPDB_FILE.exists():
    raise FileNotFoundError(
        f"CellPhoneDB database was not found at: {CPDB_FILE}"
    )


# --------------------------------------------------
# Run CellPhoneDB
# --------------------------------------------------
print("Running CellPhoneDB...")
print(f"Counts:   {COUNTS_FILE}")
print(f"Metadata: {META_FILE}")
print(f"Database: {CPDB_FILE}")
print(f"Results:  {RESULTS_DIR}")


cpdb_results = cpdb_analysis_method.call(
    cpdb_file_path=str(CPDB_FILE),
    meta_file_path=str(META_FILE),
    counts_file_path=str(COUNTS_FILE),

    # Example dataset uses human gene symbols
    counts_data="hgnc_symbol",

    # Minimum fraction of cells expressing a ligand/receptor
    threshold=0.1,

    # First test: no permutations / statistical inference
    score_interactions=False,

    output_path=str(RESULTS_DIR),

    # Easier filenames for our prototype
    output_suffix="prototype"
)


# --------------------------------------------------
# Inspect returned results
# --------------------------------------------------
print("\nCellPhoneDB completed successfully.")

print("\nReturned result objects:")
for key, value in cpdb_results.items():
    if hasattr(value, "shape"):
        print(f"  {key}: {value.shape}")
    else:
        print(f"  {key}: {type(value)}")


means = cpdb_results.get("means_result")

if means is not None:
    print("\nExample interaction rows:")
    print(
        means[
            [
                "interacting_pair",
                "gene_a",
                "gene_b",
                "receptor_a",
                "receptor_b"
            ]
        ].head(10)
    )

print("\nDone.")