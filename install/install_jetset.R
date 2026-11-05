options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("jetset", quietly = TRUE)) remotes::install_version("jetset", upgrade = "never")
cat("jetset", as.character(packageVersion("jetset")), "\n")
