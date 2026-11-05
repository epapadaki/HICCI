# Installs the R packages that conda cannot provide at the needed versions,
# into the ccc-platform environment. Safe to run again.
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

need <- function(pkg, min_version = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(TRUE)
  !is.null(min_version) && utils::packageVersion(pkg) < min_version
}

# NMF: CellChat needs >= 0.23.0; conda-forge only has 0.21.0.
# 0.28 is the version used in the environment where CellChat was validated.
if (need("NMF", "0.23.0")) {
  message("Installing NMF 0.28 from CRAN ...")
  remotes::install_version("NMF", version = "0.28", upgrade = "never", dependencies = NA)
}

# presto: fast Wilcoxon tests, used by CellChat::identifyOverExpressedGenes
if (need("presto")) {
  message("Installing presto ...")
  remotes::install_github("immunogenomics/presto", upgrade = "never", dependencies = NA)
}

# CellChat 2.2.0.9001 (jinworks/CellChat, fixed commit)
CELLCHAT_REF <- "75253cd0c9e68410e6e721a6d3a0419a1d7e358f"
if (need("CellChat")) {
  message("Installing CellChat (this compiles C++ code and can take several minutes) ...")
  remotes::install_github("jinworks/CellChat", ref = CELLCHAT_REF, upgrade = "never", dependencies = NA)
}

suppressPackageStartupMessages({ library(NMF); library(presto); library(CellChat) })
cat("R packages ready: CellChat", as.character(packageVersion("CellChat")),
    "| NMF", as.character(packageVersion("NMF")),
    "| presto", as.character(packageVersion("presto")), "\n")

# ICELLNET and its dependencies (separate script, also usable on its own)
source(file.path("install", "install_icellnet.R"))
