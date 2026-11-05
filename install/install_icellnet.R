# Installs ICELLNET (GitHub: soumelis-lab/ICELLNET, subdir "icellnet") and its
# dependencies into the current R. Safe to run again.
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")
have <- function(p) requireNamespace(p, quietly = TRUE)

# The XML package (needed by annotate and rlist) must come from conda: compiling it
# from CRAN needs the libxml2 development files (xml2-config), which R cannot find.
if (!have("XML")) stop(paste(
  "The R package XML is missing. Install it with conda first, in this environment:",
  "  mamba install -y --freeze-installed -c conda-forge r-xml", sep = "\n"))

cran <- c("remotes", "curl", "readxl", "psych", "GGally", "gplots", "ggplot2", "RColorBrewer",
          "data.table", "gridExtra", "ggthemes", "scales", "rlist", "dplyr")
miss <- cran[!vapply(cran, have, logical(1))]
if (length(miss)) { message("Installing from CRAN: ", paste(miss, collapse = ", ")); install.packages(miss) }

bioc <- c("BiocGenerics", "org.Hs.eg.db", "hgu133plus2.db", "annotate")
miss <- bioc[!vapply(bioc, have, logical(1))]
if (length(miss)) {
  if (!have("BiocManager")) install.packages("BiocManager")
  message("Installing from Bioconductor: ", paste(miss, collapse = ", "))
  BiocManager::install(miss, update = FALSE, ask = FALSE)
}

# jetset was archived on CRAN ("not available for this version of R"): install its last archived release.
if (!have("jetset")) {
  message("Installing jetset from the CRAN archive ...")
  remotes::install_version("jetset", upgrade = "never")
}

if (!have("icellnet")) {
  message("Installing ICELLNET from GitHub ...")
  remotes::install_github("soumelis-lab/ICELLNET", ref = "master", subdir = "icellnet", upgrade = "never", dependencies = NA)
}
suppressPackageStartupMessages(library(icellnet))
cat("ICELLNET ready:", as.character(packageVersion("icellnet")), "\n")
