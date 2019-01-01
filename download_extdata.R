library("dplyr")
library("RCurl")
library("XML")
library("rvest")

extdata_dir <- "inst/extdata/"
hla_alignmnets_url <- "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/"
hla_alignmnets_db <- getURL(hla_alignmnets_url,
                            verbose=TRUE,
                            ftp.use.epsv= FALSE,
                            dirlistonly = TRUE
)
hla_alignmnets_db <- strsplit(hla_alignmnets_db, "\n")[[1]]
for (aln in hla_alignmnets_db) {
  if (grepl("_prot.txt$", aln)) {
    download.file(paste0(hla_alignmnets_url, aln), destfile = paste0(extdata_dir, aln))
  }
}

aa_code_url <- "http://www.fao.org/docrep/004/y2775e/y2775e0e.htm"
aa_code_table <- html(aa_code_url) %>%
  html_table(fill=TRUE)
write.table(aa_code_table[[1]], file = "inst/extdata/aa_code.tsv", sep = "\t", col.names = FALSE, row.names = FALSE)
