#!/usr/bin/env R
# Download KIR frequencies from http://www.allelefrequencies.net/
devtools::load_all()
library("rvest")
library("xml2")

kir_gene_vector <-
  c("3DL3", "2DS2", "2DL2", "2DL3", "2DP1", "2DL1", "3DP1", "2DL4", "3DL1", "3DS1", "2DL5", "2DS3", "2DS5", "2DS4", "2DS1", "3DL2")
base_url <- "http://www.allelefrequencies.net/kir6002a.asp?"
qloci <- "kir_locus="
tab_list <- list()
for (loci in kir_gene_vector) {
  i <- 1
  query <-
    paste0(base_url, "page=", as.character(i), "&", qloci, loci)
  res <- xml2::read_html(query)
  node <-
    rvest::html_node(res, xpath = "//*[@id='divGenDetail']/table")
  if (is.na(node))
    next()
  tab <-  rvest::html_table(node)
  colnames(tab) <-
    gsub("[[:blank:]]+", " ", gsub("[^a-zA-Z%() ]", "", colnames(tab)))
  while (!is.na(node) && nrow(tab) > 0) {
    cat(loci, i, "\n")
    tab_list[[length(tab_list) + 1]] <- tab
    i <- i + 1
    query <-
      paste0(base_url, "page=", as.character(i), "&", qloci, loci)
    res <- xml2::read_html(query)
    node <-
      rvest::html_node(res, xpath = "//*[@id='divGenDetail']/table")
    tab <- rvest::html_table(node)
    colnames(tab) <-
      gsub("[[:blank:]]+", " ", gsub("[^a-zA-Z%() ]", "", colnames(tab)))
  }
}

# drop NA columns
# some columns are duplicated, and sometimes both duplicates are empty we want
# to keep one
tab_list <- lapply(
  X = tab_list,
  FUN = function(x) {
    cols <- setNames(1:ncol(x), colnames(x))
    sel <- c()
    for (i in cols) {
      is_na <- all(is.na(x[, i, drop = TRUE]))
      if (!is_na) {
        sel <- c(sel, i)
      } else {
        dup <- sum(names(cols) == names(cols)[i]) > 1
        if (dup) {
          names(cols)[i] <- ""
        } else {
          sel <- c(sel, i)
        }
      }
    }
    x[, sel]
  }
)

tab_total <- do.call(rbind, tab_list)
i_na <-
  vapply(tab_total, function(x)
    length(unique(x)) == 1, logical(1))
kirfrequency_df <- tab_total[,!i_na]
# 'Line' column specifies record number in search results
kirfrequency_df <-
  kirfrequency_df[, colnames(kirfrequency_df) != "Line"]
# here we are interested only in KIR genes, but some are alleles we remove them
sel <- grep("\\*.*", kirfrequency_df$Allele, invert = TRUE)
kirfrequency_df <- kirfrequency_df[sel,]
colnames(kirfrequency_df) <-
  c("var",
    "population",
    "proc.ind.with.allele",
    "frequency",
    "sample.size")


# order columns
kir_frequencies <-
  kirfrequency_df[, c(
    "var",
    "population",
    "proc.ind.with.allele" # here we take carrier frequency, as it is only relevant to current KIR typing implementations
  )]
colnames(kir_frequencies)[colnames(kir_frequencies) == "proc.ind.with.allele"] <- "frequency"

# set frequencies to decimals
kir_frequencies$frequency <- kir_frequencies$frequency / 100

# add KIR prefix to gene names
kir_frequencies$var <- paste0("KIR", kir_frequencies$var)

usethis::use_data(kir_frequencies, overwrite = TRUE)
