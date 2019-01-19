library("dplyr")
library("RCurl")
library("XML")
library("rvest")
devtools::load_all()

extdata_dir <- "inst/extdata/"
hla_alignmnets_url <- "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/"
hla_alignmnets_db <- getURL(hla_alignmnets_url,
                            verbose=TRUE,
                            ftp.use.epsv= FALSE,
                            dirlistonly = TRUE
)
hla_alignmnets_db <- strsplit(hla_alignmnets_db, "\n")[[1]] %>%
  grep("_prot.txt", ., value = TRUE) %>%
  grep("ClassI_prot.txt", ., invert = TRUE, value = TRUE)
for (aln in hla_alignmnets_db) {
    download.file(paste0(hla_alignmnets_url, aln), destfile = paste0(extdata_dir, aln))
}
# check if alignemnt files contains mixed records
multi_aln <- lapply(
  paste0(extdata_dir, hla_alignmnets_db),
  function(aln_file) {
    numbers <- rownames(readHlaAlignments(aln_file))
    numbers <- vapply(strsplit(numbers, "\\*"), `[[`, 1, FUN.VALUE = character(length = 1))
    return(length(unique(numbers)) > 1)
  }
) %>%
  unlist()
for (aln_file in paste0(extdata_dir, hla_alignmnets_db)[multi_aln]) {
  aln <-readHlaAlignments(aln_file)
  genes <- rownames(aln)
  ref_name <- genes[1]
  genes <- vapply(strsplit(genes, "\\*"), `[[`, 1, FUN.VALUE = character(length = 1))
  genes <- unique(genes)
  raw_aln <- readLines(aln_file)
  ref_seq_ids <- grep(ref_name, raw_aln, fixed = TRUE)
  for (i in 1:length(genes)) {
    keep <- ! grepl(paste(genes[-i], collapse = "|"), raw_aln)
    keep[ref_seq_ids] <- TRUE
    writeLines(text = raw_aln[keep],
               con = paste0(extdata_dir, genes[i], "_prot.txt"),
               sep = "\n"
    )
  }
  unlink(aln_file)
}

aprot_fasta <- "ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/A_prot.fasta"
download.file(aprot_fasta, destfile = paste0(extdata_dir, "A_prot.fasta"))
