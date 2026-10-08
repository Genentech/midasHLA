#!/usr/bin/env R
# By Migdal 2018
# Downloads HLA protein alignment files from the IPD-IMGT/HLA GitHub repository
# (https://github.com/ANHIG/IMGTHLA) for all genes shipped with the package.
# The script should be run from the package repository root, files are saved
# to the 'alignments' directory, which is then used by 'parse_alignments.R'.

# IPD-IMGT/HLA release to download, given as the first command line argument;
# the repository keeps a branch for each release, eg. "3650" for release
# 3.65.0, "Latest" (default) points to the newest one
release <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(release)) release <- "Latest"
out_dir <- "alignments"
options(timeout = max(600, getOption("timeout")))

dir.create(out_dir, showWarnings = FALSE)
genes <- sub(
  pattern = "_prot.Rdata$",
  replacement = "",
  x = list.files(file.path("inst", "extdata"), pattern = "_prot.Rdata$")
)
stopifnot(
  "no alignments found in inst/extdata, run the script from the repository root" = length(genes) > 0
)

url <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/%s/alignments/%s_prot.txt"
for (gene in genes) {
  download.file(
    url = sprintf(url, release, gene),
    destfile = file.path(out_dir, paste0(gene, "_prot.txt"))
  )
}
