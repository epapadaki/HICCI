# Run SingleCellSignalR (v1, LRscore) on the prepared expression data.
#
# Input  (prepared for CellChat, shared):  data/cellchat/
#   endometrium_normalized_expression.mtx   genes x cells, normalised
#   endometrium_genes.tsv                   column "gene"
#   endometrium_cells.tsv                   column "cell_id"
#   endometrium_metadata.tsv                columns "cell_id", "labels"
# Output: results/singlecellsignalr/
#   singlecellsignalr_raw.csv               sender, receiver, ligand, receptor, interaction type, LRscore
#   run_info.csv                            versions and parameters
#
# SingleCellSignalR v1 reports ligand-receptor pairs whose LRscore is at or
# above s.score (0.5 by default, the package's recommended threshold).
# LRscore = sqrt(l * r) / (mu + sqrt(l * r)), with l and r the mean
# expression of ligand and receptor in the two populations and mu the mean
# of the expression matrix (Cabello-Aguilar et al., Nucleic Acids Res 2020).

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellSignalR)
})

S_SCORE <- 0.5
project_root <- normalizePath(".")
input_dir <- file.path(project_root, "data", "cellchat")
output_dir <- file.path(project_root, "results", "singlecellsignalr")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("cell_signaling", envir = asNamespace("SingleCellSignalR"))) {
  stop(paste0("This script needs SingleCellSignalR version 1 (function cell_signaling). ",
              "Installed version: ", as.character(packageVersion("SingleCellSignalR")),
              ". Version 2 requires R >= 4.5 and a different interface."))
}

cat("\n=== LOADING INPUT ===\n")
m <- Matrix::readMM(file.path(input_dir, "endometrium_normalized_expression.mtx"))
genes <- read.delim(file.path(input_dir, "endometrium_genes.tsv"), stringsAsFactors = FALSE)[[1]]
cells <- read.delim(file.path(input_dir, "endometrium_cells.tsv"), stringsAsFactors = FALSE)[[1]]
meta <- read.delim(file.path(input_dir, "endometrium_metadata.tsv"), stringsAsFactors = FALSE)
if (nrow(m) != length(genes) || ncol(m) != length(cells)) stop("Matrix dimensions do not match genes and cells.")
labels <- meta$labels[match(cells, meta$cell_id)]
if (any(is.na(labels))) stop("Some cells have no cell-type label.")

# Populations as 1..k; neutral names so that population names with "-" cannot break parsing.
pops <- sort(unique(labels))
cluster <- match(labels, pops)
tokens <- paste0("P", seq_along(pops))
cat("Cells:", ncol(m), " Genes:", nrow(m), " Populations:", length(pops), "\n")

data <- as.data.frame(as.matrix(m))
rownames(data) <- genes

collect <- function(res, type) {
  if (is.null(res) || length(res) == 0) return(NULL)
  out <- list()
  for (nm in names(res)) {
    tab <- res[[nm]]
    if (is.null(tab) || nrow(tab) == 0) next
    parts <- strsplit(nm, "-", fixed = TRUE)[[1]]
    s <- pops[match(parts[1], tokens)]
    r <- pops[match(parts[length(parts)], tokens)]
    out[[length(out) + 1]] <- data.frame(
      sender_cell_type = s, receiver_cell_type = r,
      ligand = as.character(tab[[1]]), receptor = as.character(tab[[2]]),
      interaction_type = as.character(tab[["interaction type"]]),
      LRscore = as.numeric(tab[["LRscore"]]),
      mode = type, stringsAsFactors = FALSE)
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

cat("\n=== PARACRINE SIGNALLING ===\n")
para <- cell_signaling(data = data, genes = genes, cluster = cluster, int.type = "paracrine",
                       c.names = tokens, s.score = S_SCORE, species = "homo sapiens", write = FALSE, verbose = FALSE)
cat("\n=== AUTOCRINE SIGNALLING ===\n")
auto <- cell_signaling(data = data, genes = genes, cluster = cluster, int.type = "autocrine",
                       c.names = tokens, s.score = S_SCORE, species = "homo sapiens", write = FALSE, verbose = FALSE)

res <- rbind(collect(para, "paracrine"), collect(auto, "autocrine"))
if (is.null(res) || nrow(res) == 0) stop("SingleCellSignalR returned no interactions.")
res <- res[!duplicated(res[, c("sender_cell_type", "receiver_cell_type", "ligand", "receptor")]), ]
write.csv(res, file.path(output_dir, "singlecellsignalr_raw.csv"), row.names = FALSE)

info <- data.frame(
  field = c("method", "package_version", "r_version", "s_score", "species", "modes", "n_cells", "n_genes",
            "n_populations", "n_interactions", "input"),
  value = c("SingleCellSignalR", as.character(packageVersion("SingleCellSignalR")),
            paste(R.version$major, R.version$minor, sep = "."), S_SCORE, "homo sapiens",
            "paracrine; autocrine", ncol(m), nrow(m), length(pops), nrow(res),
            "normalised expression prepared for CellChat"))
write.csv(info, file.path(output_dir, "run_info.csv"), row.names = FALSE)
cat("\nSingleCellSignalR finished:", nrow(res), "interactions (LRscore >=", S_SCORE, ")\n")
