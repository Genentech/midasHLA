context("parsing functions")

test_that("readHlaCalls", {
  file <- system.file("extdata", "MiDAS_tut_HLA.txt", package = "midasHLA")
  hla_calls <- readHlaCalls(file)
  load(system.file("extdata", "test_hla_calls.Rdata", package = "midasHLA"))
  expect_equal(hla_calls, test_hla_calls)

  hla_calls_res2 <- readHlaCalls(file, resolution = 2)
  res2 <- getAlleleResolution(unlist(hla_calls_res2[, -1]))
  load(system.file("extdata", "test_hla_calls_res.Rdata", package = "midasHLA"))
  expect_equal(res2, test_hla_calls_res)

  expect_error(readHlaCalls(file.path("path", "to", "nonexisting", "file")),
               sprintf("Path '%s' does not exist",
                       file.path("path", "to", "nonexisting", "file")
               )
  )

  expect_error(readHlaCalls(file, resolution = "foo"),
               "resolution is not a count \\(a single positive integer\\)"
  )

  expect_error(readHlaCalls(file, na.strings = 1),
               "na.strings is not a character vector"
  )

  fake_calls <- data.frame(ID = c("Sample1", "Sample2", "Sample3"),
                           A_1 = c("A*01", "A*02", "A*03"),
                           A_2 = c("A*01", "B*02", "C*03")
  )

  fake_calls_non_uniq_genes <- tempfile()
  write.table(fake_calls,
              file = fake_calls_non_uniq_genes,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE
  )
  expect_error(readHlaCalls(fake_calls_non_uniq_genes, resolution = 2),
               "Gene names in columns are not identical"
  )
  unlink(fake_calls_non_uniq_genes)

  fake_calls_NA_col <- tempfile()
  fake_calls[, 3] <- NA
  write.table(fake_calls,
              file = fake_calls_NA_col,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE
  )
  expect_error(readHlaCalls(fake_calls_NA_col, resolution = 2),
               "One of the columns contains only NA")
  unlink(fake_calls_NA_col)
})

test_that("readHlaAlignments handles gaps in reference allele", {
  # since IPD-IMGT/HLA 3.5x the reference allele contains gaps ('.') in columns
  # where other alleles have insertions; these columns are not numbered, they
  # are named after the preceding numbered position, eg. 4.1, 4.2
  file <- test_path("test_ref_gaps_prot.txt")
  aln <- readHlaAlignments(file)
  pos <- c("-3", "-2", "-2.1", "-1", "1", "1.1", "2", "3", "4", "4.1", "4.2",
           "5", "6", "7", "8", "9", "10", "11", "11.1", "12", "13")
  ref <- c("M", "K", ".", "A", "G", ".", "S", "H", "S", ".", ".", "M", "R",
           "Y", "F", "F", "T", "S", ".", "V", "S")

  expect_equal(colnames(aln), pos)
  expect_equal(
    rownames(aln),
    c("X*01:01:01:01", "X*01:02", "X*02:01", "X*03:01", "X*04:01")
  )
  expect_equal(unname(aln["X*01:01:01:01", ]), ref)
  # insertions
  expect_equal(unname(aln["X*01:02", ]), replace(ref, pos == "1.1", "K"))
  expect_equal(
    unname(aln["X*02:01", ]),
    replace(ref, pos %in% c("4.1", "4.2", "9"), c("R", "R", "Y"))
  )
  expect_equal(
    unname(aln["X*04:01", ]),
    replace(ref, pos %in% c("-2.1", "13"), c("Q", "."))
  )
  # deletion
  expect_equal(unname(aln["X*03:01", ]), replace(ref, pos == "3", "."))

  aln_trim <- readHlaAlignments(file, trim = TRUE)
  expect_equal(aln_trim, aln[, which(pos == "1"):length(pos)])
})

test_that("readHlaAlignments agrees with IPD-IMGT/HLA protein records", {
  # alignments/<gene>_prot.txt are excerpts of IPD-IMGT/HLA 3.65.0 alignment
  # files; expected alignments are built from IPD-IMGT/HLA REST API records
  # (see helper-ipd.R), which are independent of the alignment files
  for (gene in c("A", "B", "C", "DQB1")) {
    ipd <- readIpdProtein(gene)
    aln <- readHlaAlignments(
      test_path("alignments", paste0(gene, "_prot.txt")),
      unkchar = "*"
    )
    expected <- ipdAlignment(ipd)

    # columns beyond allele sequence: no residue
    aln[aln == ""] <- "."
    expected[expected == ""] <- "."

    expect_equal(colnames(aln), colnames(expected), info = gene)
    expect_equal(aln[ipd$allele, ], expected, info = gene)

    for (i in seq_len(nrow(ipd))) {
      allele <- ipd$allele[i]
      # parsing is lossless: residues read in order give allele's sequence
      residues <- aln[allele, ! aln[allele, ] %in% c(".", "*")]
      expect_equal(
        paste(residues, collapse = ""),
        ipd$protein[i],
        info = allele
      )
      # mature protein of allele starts at position 1
      if (! is.na(ipd$signal_length[i])) {
        signal <- as.integer(ipd$signal_length[i])
        expect_equal(
          aln[allele, "1"],
          substr(ipd$protein[i], signal + 1, signal + 1),
          info = allele
        )
      }
    }
  }

  # known insertions
  aln <- readHlaAlignments(test_path("alignments", "C_prot.txt"))
  expect_equal(
    unname(aln["C*17:03:01:01", paste0("300.", 1:6)]),
    c("A", "V", "L", "A", "V", "L")
  )
  expect_equal(
    unname(aln["C*01:02:01:01", paste0("300.", 1:6)]),
    rep(".", 6)
  )
  aln <- readHlaAlignments(test_path("alignments", "B_prot.txt"))
  expect_equal(unname(aln["B*73:01:01:01", "296.1"]), "T")
})

test_that("readHlaAlignments", {
  file <- system.file("extdata", "TAP1_prot.txt", package = "midasHLA")
  # our alignment files also contain infered alignment sequences
  m <- c("TAP1*01:02N", "TAP1*03:01", "TAP1*04:01", "TAP1*05:01", "TAP1*06:01", 
         "TAP1*02:01:01", "TAP1*02:01:02", "TAP1*01:01:01:01", "TAP1*01:01:01:02", 
         "TAP1*01:01:01:03", "TAP1*01:01:01:04", "TAP1*01:01:01:05")
  hla_alignments <- readHlaAlignments(file)[m, ]
  test_hla_alignments <- readHlaAlignments(gene = "TAP1")[m, ]
  expect_equal(hla_alignments, test_hla_alignments)

  expect_error(readHlaAlignments(file.path("path", "to", "nonexisting", "file")),
               sprintf("Path '%s' does not exist",
                       file.path("path", "to", "nonexisting", "file")
               )
  )

  expect_error(
    readHlaAlignments(system.file("extdata", "A_prot.txt", package = "midasHLA"),
                      trim = c(TRUE, TRUE),
                      unkchar = ""
    ),
    "trim is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    readHlaAlignments(system.file("extdata", "A_prot.txt", package = "midasHLA"),
                      trim = TRUE,
                      unkchar = c("a", "b", "c")
    ),
    "unkchar is not a string \\(a length one character vector\\)."
  )

  expect_error(readHlaAlignments(gene = "foo"),
               "alignment for FOO is not available"
  )

  aln_file <- system.file("extdata/TAP1_prot.txt", package = "midasHLA")
  hla_alignments <- readHlaAlignments(aln_file, trim = FALSE)
  fasta_file <- system.file("extdata", "TAP1_prot.fasta", package = "midasHLA")
  fasta <- seqinr::read.alignment(fasta_file, format = "fasta")
  fasta <- fasta$seq[[1]]
  expect_equal(paste(hla_alignments["TAP1*01:01:01:01", ], collapse = ""),
               toupper(fasta)
  ) # check sequence with sequence from fasta

  hla_alignments_trim <- readHlaAlignments(aln_file, trim = TRUE)
  n_trimmed <- ncol(hla_alignments) - ncol(hla_alignments_trim)
  expect_equal(n_trimmed, 0) # check if sequence is trimmed properly

  fake_aln <- readLines(aln_file)
  fake_aln <- vapply(X = fake_aln,
                     FUN = function(x) {
                       number <- stri_split_regex(x, "\\s+")[[1]]
                       ifelse(any(checkAlleleFormat(number)), "empty ", x)
                     },
                     FUN.VALUE = character(length = 1),
                     USE.NAMES = FALSE
  )
  fake_aln_tmp <- tempfile()
  writeLines(text = fake_aln, con = fake_aln_tmp)
  expect_error(readHlaAlignments(fake_aln_tmp),
               "could not find alleles numbers in the alignment file"
  )
  unlink(fake_aln_tmp)

  # alignments lines contain non standard characters
  fake_aln <- readLines(aln_file)
  fake_aln <- vapply(X = fake_aln,
                     FUN = function(x) {
                       number <- stri_split_regex(x, "\\s+")[[1]]
                       if (any(checkAlleleFormat(number))) {
                         li <- length(number) - 1
                         number[li] <- paste(sample(c("?", ">", "<", "#", "@", "!"),
                                                    nchar(number[li]),
                                                    replace = TRUE
                                             ),
                                             collapse = ""
                                       )
                         paste(number, collapse = " ")
                       } else {
                         x
                       }
                     },
                     FUN.VALUE = character(length = 1),
                     USE.NAMES = FALSE
  )
  fake_aln_tmp <- tempfile()
  writeLines(text = fake_aln, con = fake_aln_tmp)
  expect_error(readHlaAlignments(fake_aln_tmp),
               "alignments lines contain non standard characters"
  )
  unlink(fake_aln_tmp)

  fake_aln <- readLines(aln_file)
  fake_aln[grepl("Prot", fake_aln)] <- "Prot"
  fake_aln_tmp <- tempfile()
  writeLines(text = fake_aln, con = fake_aln_tmp)
  expect_error(readHlaAlignments(fake_aln_tmp, trim=TRUE),
               "start codon is not marked properly in the input file"
  )
  unlink(fake_aln_tmp)
})

test_that("readKirCalls", {
  kpi_output <- data.frame(
    ID = c("SAM24320917", "SAM24320918"),
    KIR3DL3 = c(1L, 1L),
    KIR2DS2 = 0:1,
    KIR2DL2 = 0:1,
    KIR2DL3 = c(1L, 1L),
    KIR2DP1 = c(1L, 1L),
    KIR2DL1 = c(1L, 1L),
    KIR3DP1 = c(1L, 1L),
    KIR2DL4 = c(1L, 1L),
    KIR3DL1 = c(1L, 1L),
    KIR3DS1 = c(0L, 0L),
    KIR2DL5 = c(0L, 0L),
    KIR2DS3 = c(0L, 0L),
    KIR2DS5 = c(0L, 0L),
    KIR2DS4 = c(1L, 1L),
    KIR2DS1 = 1:0,
    KIR3DL2 = c(1L, 1L),
    stringsAsFactors = FALSE
  )
  file <- tempfile()
  write.table(
    x = kpi_output,
    file = file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  kir_calls <- readKirCalls(file)
  test_kir_calls <- kpi_output
  expect_equal(kir_calls, test_kir_calls)

  colnames(kpi_output)[1] <- "SAMID"
  write.table(
    x = kpi_output,
    file = file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  expect_error(readKirCalls(file), "Columns: 'SAMID' in kir_calls should be named 'ID'")

  kpi_output <- kpi_output[-1]
  write.table(
    x = kpi_output,
    file = file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
  expect_error(readKirCalls(file), "kir_calls shiuld have 17 columns: ID, KIR3DL3, KIR2DS2, KIR2DL2, KIR2DL3, KIR2DP1, KIR2DL1, KIR3DP1, KIR2DL4, KIR3DL1, KIR3DS1, KIR2DL5, KIR2DS3, KIR2DS5, KIR2DS4, KIR2DS1, KIR3DL2")
  unlink(file)
})
