#!/usr/bin/env R
# Download HLA frequencies from http://www.allelefrequencies.net/
devtools::load_all()
library("rvest")
library("xml2")

midas_var_hla_locus_vector <-
  c("A", "B", "C", "DQA1", "DQB1", "DPA1", "DPB1", "DRB1", "DRB3", "DRB4", "DRB5")
if (! file.exists("inst/extdata/allelefrequency_df.Rdata")) {
  base_url <- "http://www.allelefrequencies.net/hla6006a.asp?"
  qloci <- "hla_locus="
  tab_list <- list()
  for (loci in midas_var_hla_locus_vector) {
    i <- 1
    query <- paste0(base_url, "page=", as.character(i), "&", qloci, loci)
    res <- xml2::read_html(query)
    node <- rvest::html_node(res, xpath = "//*[@id='divGenDetail']/table")
    if (is.na(node)) next()
    tab <-  rvest::html_table(node)
    colnames(tab) <-
      gsub("[[:blank:]]+", " ", gsub("[^a-zA-Z%() ]", "", colnames(tab)))
    while (! is.na(node) && nrow(tab) > 0) {
      tab_list[[length(tab_list) + 1]] <- tab
      i <- i + 1
      query <- paste0(base_url, "page=", as.character(i), "&", qloci, loci)
      res <- xml2::read_html(query)
      node <- rvest::html_node(res, xpath = "//*[@id='divGenDetail']/table")
      tab <- rvest::html_table(node)
      colnames(tab) <-
        gsub("[[:blank:]]+", " ", gsub("[^a-zA-Z%() ]", "", colnames(tab)))
    }
  }

  tab_total <- do.call(rbind, tab_list)
  # Some cols contains NA other links to external website
  i_na <- vapply(tab_total, function(x) length(unique(x)) == 1, logical(1))
  allelefrequency_df <- tab_total[, ! i_na]
  # 'Line' column specifies record number in search results
  allelefrequency_df <-
    allelefrequency_df[, colnames(allelefrequency_df) != "Line"]
  colnames(allelefrequency_df) <-
    c("var", "population", "proc.ind.with.allele", "frequency", "sample.size")

  # Allele Frequencies marked with (*) were calculated from all alleles in the corresponding G group.
  # These symbols can be found in "proc.ind.with.allele", "frequency" columns
  star_cols <- c("proc.ind.with.allele", "frequency")
  new_star_cols <-
    c(
      "proc.ind.with.allele",
      "proc.ind.with.allele.ggroup",
      "frequency",
      "frequency.ggroup"
    )
  ggroup_cols <- c("proc.ind.with.allele.ggroup", "frequency.ggroup")
  data.table::setDT(allelefrequency_df)
  allelefrequency_df[,
                     (new_star_cols) := unlist(
                       lapply(
                         X = .SD,
                         FUN = function(x) {
                           data.table::tstrsplit(x,
                                                 split = " ",
                                                 fixed = TRUE,
                                                 type.convert = TRUE
                           )
                         }
                       ),
                       recursive = FALSE
                     ),
                     .SDcols =  star_cols
                     ][,
                       (ggroup_cols) := lapply(
                         X = .SD,
                         FUN = function(x) {
                           ! is.na(x)
                         }
                       ),
                       .SDcols = ggroup_cols
                       ][,
                         sample.size := as.numeric(gsub(",", "", sample.size))
                         ]
  data.table::setDF(allelefrequency_df)
  save(allelefrequency_df, file = "inst/extdata/allelefrequency_df.Rdata") # this object is save for development purposes to avoid scrapping the database over and over
} else {
  load(system.file("extdata", "allelefrequency_df.Rdata", package = "MiDAS"))
}

# order columns
allele_frequencies <-
  allelefrequency_df[, c(
    "var",
    "population",
    "frequency"
    # "frequency.ggroup",
    # "sample.size",
    # "proc.ind.with.allele",
    # "proc.ind.with.allele.ggroup"
  )]

usethis::use_data(allele_frequencies, overwrite = TRUE)
