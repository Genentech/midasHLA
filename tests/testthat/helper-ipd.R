# Helpers building expected protein alignments from IPD-IMGT/HLA REST API
# records (https://www.ebi.ac.uk/cgi-bin/ipd/api/allele/<accession>), stored in
# alignments/ipd_protein.tsv. Each record holds allele's protein sequence,
# signal peptide length and protein CIGAR string placing allele's residues in
# the IPD-IMGT/HLA alignment (M - residue, D - no residue, N - not sequenced).
# The data is independent of the alignment text files read by
# readHlaAlignments, so it can be used to validate them.

readIpdProtein <- function(gene) {
  ipd <- read.table(
    file = test_path("alignments", "ipd_protein.tsv"),
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  ipd[ipd$gene == gene, ]
}

expandProteinCigar <- function(cigar, protein) {
  ops <- regmatches(cigar, gregexpr("[0-9]*[MDN]", cigar))[[1]]
  assertthat::assert_that(
    identical(paste(ops, collapse = ""), cigar),
    msg = "unexpected CIGAR operation"
  )
  residues <- strsplit(protein, "")[[1]]
  out <- character(0)
  i <- 0
  for (op in ops) {
    n <- suppressWarnings(as.integer(substr(op, 1, nchar(op) - 1)))
    if (is.na(n)) n <- 1
    out <- c(out, switch(
      substr(op, nchar(op), nchar(op)),
      M = residues[i + seq_len(n)],
      D = rep(".", n),
      N = rep("*", n)
    ))
    if (endsWith(op, "M")) i <- i + n # only residues are present in protein
  }

  out
}

# Expected alignment of alleles in ipd, with first row being the reference
# allele. Columns are named according to IPD-IMGT/HLA numbering of the
# reference allele: residues are numbered from -signal_length, skipping 0;
# columns without residue in the reference (insertions) are named after the
# preceding numbered position, eg. 4.1, 4.2. Columns beyond allele sequence
# are marked as "".
ipdAlignment <- function(ipd) {
  cols <- Map(expandProteinCigar, ipd$protein_cigar, ipd$protein)
  width <- max(lengths(cols))
  aln <- do.call(rbind, lapply(cols, function(x) c(x, rep("", width - length(x)))))
  rownames(aln) <- ipd$allele

  ref <- aln[1, ]
  pos <- character(width)
  n <- -as.integer(ipd$signal_length[1]) - 1
  k <- 0
  for (i in seq_len(width)) {
    if (ref[i] %in% c(".", "")) {
      k <- k + 1
      pos[i] <- paste0(n, ".", k)
    } else {
      n <- if (n == -1) 1 else n + 1
      k <- 0
      pos[i] <- as.character(n)
    }
  }
  colnames(aln) <- pos

  aln
}
