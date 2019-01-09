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

aprot_fasta <- "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/A_prot.fasta"
download.file(aprot_fasta, destfile = paste0(extdata_dir, "A_prot.fasta"))
