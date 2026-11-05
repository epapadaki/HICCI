# Run ICELLNET (v2) on the prepared expression data, for every sender -> receiver pair.
#
# Input  (prepared for CellChat, shared): data/cellchat/  (normalised matrix, genes, cells, labels)
#        ICELLNET database: data/icellnet/ICELLNETdb.tsv  (local snapshot; downloaded if absent)
# Output: results/icellnet/icellnet_raw.csv, results/icellnet/run_info.csv
#
# Following the ICELLNET scRNA-seq tutorial (Example 2): average normalised
# expression per population, genes expressed in fewer than FILTER_PERC % of a
# population's cells set to 0 (sc.data.cleaning behaviour), gene.scaling()
# (maximum of each gene = 10), then icellnet.score(direction = "out") with each
# population as central cell. The score of an interaction is the product of
# the scaled ligand and receptor expression (0-100); complexes are handled by
# ICELLNET's own rules. Interactions with a score of 0 are not reported.

suppressPackageStartupMessages({
  library(Matrix)
  library(icellnet)
})

FILTER_PERC <- 10
DB_URL <- "https://raw.githubusercontent.com/soumelis-lab/ICELLNET/master/data/ICELLNETdb.tsv"
project_root <- normalizePath(".")
input_dir <- file.path(project_root, "data", "cellchat")
db_file <- file.path(project_root, "data", "icellnet", "ICELLNETdb.tsv")
output_dir <- file.path(project_root, "results", "icellnet")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(db_file)) {
  cat("Downloading the ICELLNET database ...\n")
  dir.create(dirname(db_file), recursive = TRUE, showWarnings = FALSE)
  utils::download.file(DB_URL, db_file, quiet = TRUE)
}
db <- as.data.frame(read.csv(db_file, sep = "\t", header = TRUE, check.names = FALSE,
                             stringsAsFactors = FALSE, na.strings = ""))
cat("ICELLNET database:", nrow(db), "interactions\n")

cat("\n=== LOADING INPUT ===\n")
m <- Matrix::readMM(file.path(input_dir, "endometrium_normalized_expression.mtx"))
m <- as(m, "CsparseMatrix")
genes <- read.delim(file.path(input_dir, "endometrium_genes.tsv"), stringsAsFactors = FALSE)[[1]]
cells <- read.delim(file.path(input_dir, "endometrium_cells.tsv"), stringsAsFactors = FALSE)[[1]]
meta <- read.delim(file.path(input_dir, "endometrium_metadata.tsv"), stringsAsFactors = FALSE)
labels <- meta$labels[match(cells, meta$cell_id)]
if (any(is.na(labels))) stop("Some cells have no cell-type label.")
pops <- sort(unique(labels))
tokens <- paste0("P", seq_along(pops))       # neutral names: ICELLNET works on data.frame column names
cat("Cells:", ncol(m), " Genes:", nrow(m), " Populations:", length(pops), "\n")

# Average expression and share of expressing cells per population
ind <- Matrix::sparseMatrix(i = seq_along(labels), j = match(labels, pops), x = 1,
                            dims = c(length(labels), length(pops)))
n <- Matrix::colSums(ind)
avg <- as.matrix(m %*% ind) %*% diag(1 / n)
pct <- as.matrix((m > 0) %*% ind) %*% diag(1 / n)
if (FILTER_PERC > 0) avg[pct < FILTER_PERC / 100] <- 0
avg <- as.data.frame(avg)
colnames(avg) <- tokens
rownames(avg) <- genes

data.icell <- as.data.frame(gene.scaling(avg, n = 1, db = db))
if (!"Symbol" %in% colnames(data.icell)) data.icell$Symbol <- rownames(data.icell)

cat("\n=== SCORING (each population as central cell) ===\n")
PC.data <- as.data.frame(data.icell[, c(tokens, "Symbol")], row.names = rownames(data.icell))
out <- list()
for (k in seq_along(tokens)) {
  cc <- tokens[k]
  res <- icellnet.score(direction = "out", PC.data = PC.data,
                        CC.data = as.data.frame(data.icell[, cc], row.names = rownames(data.icell)),
                        PC = tokens, CC.type = "RNAseq", PC.type = "RNAseq", db = db)
  lr <- as.matrix(res[[2]])
  lr[is.na(lr)] <- 0
  idx <- which(lr > 0, arr.ind = TRUE)
  if (nrow(idx) == 0) next
  out[[length(out) + 1]] <- data.frame(
    sender_cell_type = pops[k],
    receiver_cell_type = pops[match(colnames(lr)[idx[, 2]], tokens)],
    pair_name = rownames(lr)[idx[, 1]],
    score = lr[idx],
    stringsAsFactors = FALSE)
  cat(" ", pops[k], ":", nrow(idx), "interactions\n")
}
res <- do.call(rbind, out)
if (is.null(res) || nrow(res) == 0) stop("ICELLNET returned no interactions.")
write.csv(res, file.path(output_dir, "icellnet_raw.csv"), row.names = FALSE)

info <- data.frame(
  field = c("method", "package_version", "r_version", "filter_perc", "scaling", "direction", "database",
            "database_interactions", "n_cells", "n_genes", "n_populations", "n_interactions"),
  value = c("ICELLNET", as.character(packageVersion("icellnet")), paste(R.version$major, R.version$minor, sep = "."),
            FILTER_PERC, "gene.scaling(n = 1): maximum per gene = 10", "out (each population as central cell)",
            "ICELLNETdb.tsv", nrow(db), ncol(m), nrow(m), length(pops), nrow(res)))
write.csv(info, file.path(output_dir, "run_info.csv"), row.names = FALSE)
cat("\nICELLNET finished:", nrow(res), "interactions with a score above 0\n")
