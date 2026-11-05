suppressPackageStartupMessages({
    library(CellChat)
    library(Matrix)
    library(future)
})


options(
    stringsAsFactors = FALSE
)


# --------------------------------------------------
# Project paths
# --------------------------------------------------
project_root <- normalizePath(
    ".",
    mustWork = TRUE
)

# Real dataset id: this script runs inside datasets/<id>/runs/<slug>/, so the
# run folder's grandparent is the dataset id (example id outside a run folder).
dataset_id <- if (basename(dirname(project_root)) == "runs") {
    basename(dirname(dirname(project_root)))
} else {
    "endometrium_example"
}


input_dir <- file.path(
    project_root,
    "data",
    "cellchat"
)


output_dir <- file.path(
    project_root,
    "results",
    "cellchat"
)


dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


matrix_file <- file.path(
    input_dir,
    "endometrium_normalized_expression.mtx"
)


genes_file <- file.path(
    input_dir,
    "endometrium_genes.tsv"
)


cells_file <- file.path(
    input_dir,
    "endometrium_cells.tsv"
)


metadata_file <- file.path(
    input_dir,
    "endometrium_metadata.tsv"
)


# --------------------------------------------------
# Validate input files
# --------------------------------------------------
required_files <- c(
    matrix_file,
    genes_file,
    cells_file,
    metadata_file
)


missing_files <- required_files[
    !file.exists(required_files)
]


if (length(missing_files) > 0) {

    stop(
        paste(
            "Missing CellChat input file(s):",
            paste(
                missing_files,
                collapse = ", "
            )
        )
    )
}


# --------------------------------------------------
# Load normalized expression
# --------------------------------------------------
cat(
    "\n=== LOADING CELLCHAT INPUT ===\n"
)


data.input <- Matrix::readMM(
    matrix_file
)


data.input <- as(
    data.input,
    "dgCMatrix"
)


genes <- read.delim(
    genes_file,
    header = TRUE,
    stringsAsFactors = FALSE
)


cells <- read.delim(
    cells_file,
    header = TRUE,
    stringsAsFactors = FALSE
)


meta <- read.delim(
    metadata_file,
    header = TRUE,
    stringsAsFactors = FALSE
)


# --------------------------------------------------
# Validate dimensions
# --------------------------------------------------
if (
    nrow(data.input)
    != nrow(genes)
) {

    stop(
        "Expression matrix rows do not match the gene table."
    )
}


if (
    ncol(data.input)
    != nrow(cells)
) {

    stop(
        "Expression matrix columns do not match the cell table."
    )
}


if (
    nrow(meta)
    != nrow(cells)
) {

    stop(
        "Metadata rows do not match the cell table."
    )
}


# --------------------------------------------------
# Assign biological names
# --------------------------------------------------
rownames(
    data.input
) <- genes$gene


colnames(
    data.input
) <- cells$cell_id


rownames(
    meta
) <- meta$cell_id


# --------------------------------------------------
# Explicit metadata alignment
# --------------------------------------------------
meta <- meta[
    colnames(data.input),
    ,
    drop = FALSE
]


if (
    any(
        is.na(meta$labels)
    )
) {

    stop(
        "Missing cell labels were found after metadata alignment."
    )
}


# --------------------------------------------------
# Expression validation
#
# IMPORTANT:
# adata.X was already normalized/transformed.
# We do NOT normalize or log-transform again.
# --------------------------------------------------
if (
    length(data.input@x) == 0
) {

    stop(
        "Expression matrix contains no non-zero values."
    )
}


if (
    min(data.input@x) < 0
) {

    stop(
        "Negative expression values detected."
    )
}


cat(
    "Expression matrix:",
    nrow(data.input),
    "genes x",
    ncol(data.input),
    "cells\n"
)


cat(
    "Expression range:",
    min(data.input@x),
    "to",
    max(data.input@x),
    "\n"
)


cat(
    "Cell types:",
    length(
        unique(meta$labels)
    ),
    "\n"
)


cat(
    "\nCell-type counts:\n"
)


print(
    table(meta$labels)
)


# --------------------------------------------------
# Create CellChat object
# --------------------------------------------------
cat(
    "\n=== CREATING CELLCHAT OBJECT ===\n"
)


cellchat <- createCellChat(
    object = data.input,
    meta = meta,
    group.by = "labels"
)


cat(
    "\nCellChat groups:\n"
)


print(
    table(cellchat@idents)
)


# --------------------------------------------------
# Configure human CellChatDB
# --------------------------------------------------
cat(
    "\n=== CONFIGURING CELLCHATDB ===\n"
)


CellChatDB.use <- subsetDB(
    CellChatDB.human
)


cellchat@DB <- CellChatDB.use


cat(
    "CellChatDB interactions used:",
    nrow(
        CellChatDB.use$interaction
    ),
    "\n"
)


cat(
    "\nDatabase annotation counts:\n"
)


print(
    table(
        CellChatDB.use$interaction$annotation
    )
)


# --------------------------------------------------
# Save exact CellChatDB snapshot used
# --------------------------------------------------
database_file <- file.path(
    output_dir,
    "cellchatdb_interactions_used.csv"
)


write.csv(
    CellChatDB.use$interaction,
    database_file,
    row.names = FALSE
)


# --------------------------------------------------
# Subset signaling genes
# --------------------------------------------------
cat(
    "\n=== SUBSETTING SIGNALING GENES ===\n"
)


cellchat <- subsetData(
    cellchat
)


cat(
    "Signaling genes retained:",
    nrow(
        cellchat@data.signaling
    ),
    "\n"
)


# --------------------------------------------------
# Parallel-processing configuration
# --------------------------------------------------
options(
    future.globals.maxSize =
        8 * 1024^3
)


future::plan(
    "multisession",
    workers = 4
)


# --------------------------------------------------
# Identify overexpressed genes
# --------------------------------------------------
cat(
    "\n=== IDENTIFYING OVEREXPRESSED GENES ===\n"
)


cellchat <- identifyOverExpressedGenes(
    cellchat
)


# --------------------------------------------------
# Identify overexpressed interactions
# --------------------------------------------------
cat(
    "\n=== IDENTIFYING OVEREXPRESSED INTERACTIONS ===\n"
)


cellchat <- identifyOverExpressedInteractions(
    cellchat
)


# --------------------------------------------------
# Infer communication probabilities
#
# triMean is intentionally conservative.
# --------------------------------------------------
cat(
    "\n=== COMPUTING COMMUNICATION PROBABILITIES ===\n"
)


set.seed(
    1
)


cellchat <- computeCommunProb(
    cellchat,
    type = "triMean",
    nboot = 100,
    seed.use = 1
)


# --------------------------------------------------
# Filter very small cell groups
# --------------------------------------------------
cat(
    "\n=== FILTERING COMMUNICATIONS ===\n"
)


cellchat <- filterCommunication(
    cellchat,
    min.cells = 10
)


# --------------------------------------------------
# Extract LR-level interactions
# --------------------------------------------------
cat(
    "\n=== EXTRACTING LR INTERACTIONS ===\n"
)


df.net <- subsetCommunication(
    cellchat,
    thresh = 0.05
)


if (
    is.null(df.net)
) {

    stop(
        "subsetCommunication returned NULL."
    )
}


cat(
    "Significant CellChat rows:",
    nrow(df.net),
    "\n"
)


cat(
    "\nOutput columns:\n"
)


print(
    colnames(df.net)
)


# --------------------------------------------------
# Save raw LR-level predictions
# --------------------------------------------------
raw_file <- file.path(
    output_dir,
    "cellchat_interactions_raw.csv"
)


write.csv(
    df.net,
    raw_file,
    row.names = FALSE
)


# --------------------------------------------------
# Compute pathway-level probabilities
# --------------------------------------------------
cat(
    "\n=== COMPUTING PATHWAY NETWORK ===\n"
)


cellchat <- computeCommunProbPathway(
    cellchat
)


# --------------------------------------------------
# Aggregate network
# --------------------------------------------------
cellchat <- aggregateNet(
    cellchat
)


# --------------------------------------------------
# Save complete CellChat object
# --------------------------------------------------
object_file <- file.path(
    output_dir,
    "cellchat_endometrium.rds"
)


saveRDS(
    cellchat,
    object_file
)


# --------------------------------------------------
# Save run provenance
# --------------------------------------------------
run_info <- data.frame(

    field = c(
        "CellChat_version",
        "dataset_id",
        "species",
        "expression_input",
        "expression_pre_normalized",
        "normalization_in_R",
        "database",
        "database_scope",
        "expression_average",
        "nboot",
        "seed",
        "pvalue_threshold",
        "min_cells",
        "cells",
        "genes",
        "cell_types",
        "significant_LR_rows"
    ),

    value = c(
        as.character(
            packageVersion("CellChat")
        ),

        dataset_id,

        "human",

        "AnnData X",

        "TRUE",

        "FALSE",

        "CellChatDB v2 human",

        paste(
            "Secreted Signaling;",
            "ECM-Receptor;",
            "Cell-Cell Contact;",
            "Non-protein excluded"
        ),

        "triMean",

        "100",

        "1",

        "0.05",

        "10",

        ncol(data.input),

        nrow(data.input),

        length(
            unique(meta$labels)
        ),

        nrow(df.net)
    ),

    stringsAsFactors = FALSE
)


run_info_file <- file.path(
    output_dir,
    "run_info.csv"
)


write.csv(
    run_info,
    run_info_file,
    row.names = FALSE
)


# --------------------------------------------------
# Restore sequential future backend
# --------------------------------------------------
future::plan(
    "sequential"
)


# --------------------------------------------------
# Final summary
# --------------------------------------------------
cat(
    "\n========================================\n"
)


cat(
    "=== CELLCHAT COMPLETE ===\n"
)


cat(
    "========================================\n"
)


cat(
    "CellChat version:",
    as.character(
        packageVersion("CellChat")
    ),
    "\n"
)


cat(
    "CellChatDB interactions used:",
    nrow(
        CellChatDB.use$interaction
    ),
    "\n"
)


cat(
    "Significant CellChat LR rows:",
    nrow(df.net),
    "\n"
)


cat(
    "\nSaved files:\n"
)


cat(
    raw_file,
    "\n"
)


cat(
    object_file,
    "\n"
)


cat(
    run_info_file,
    "\n"
)


cat(
    database_file,
    "\n"
)