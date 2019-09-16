library(assertthat)
library(testthat)

context("HLA allele information parsing")


test_that("HLA allele calls are read properly", {
  file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
  hla_calls <- readHlaCalls(file)
  load(system.file("extdata", "test_hla_calls.Rdata", package = "MiDAS"))
  expect_equal(hla_calls, test_hla_calls)

  hla_calls_res2 <- readHlaCalls(file, resolution = 2)
  res2 <- getAlleleResolution(unlist(hla_calls_res2[, -1]))
  load(system.file("extdata", "test_hla_calls_res.Rdata", package = "MiDAS"))
  expect_equal(res2, test_res2)

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

test_that("HLA allele alignments are read properly", {
  file <- system.file("extdata", "TAP1_prot.txt", package = "MiDAS")
  hla_alignments <- readHlaAlignments(file)
  load(system.file("extdata", "test_hla_alignments.Rdata", package = "MiDAS"))
  expect_equal(hla_alignments, test_hla_alignments)

  hla_alignments_res2 <- readHlaAlignments(file, resolution = 2)
  res2 <- getAlleleResolution(rownames(hla_alignments_res2))
  test_res2 <- c(2, 4, 2, 2, 2, 2, 2)
  expect_equal(res2, test_res2)

  expect_error(readHlaAlignments(file.path("path", "to", "nonexisting", "file")),
               sprintf("Path '%s' does not exist",
                       file.path("path", "to", "nonexisting", "file")
               )
  )

  expect_error(
    readHlaAlignments(system.file("extdata", "A_prot.txt", package = "MiDAS"),
                      trim = c(T, T),
                      unkchar = ""
    ),
    "trim is not a flag \\(a length one logical vector\\)."
  )

  expect_error(
    readHlaAlignments(system.file("extdata", "A_prot.txt", package = "MiDAS"),
                      trim = T,
                      unkchar = c("a", "b", "c")
    ),
    "unkchar is not a string \\(a length one character vector\\)."
  )

  expect_error(readHlaAlignments(gene = "foo"),
               "alignment for FOO is not available"
  )

  aln_file <- system.file("extdata/A_prot.txt", package = "MiDAS")
  hla_alignments <- readHlaAlignments(aln_file, trim = FALSE)
  fasta_file <- system.file("extdata", "A_prot.fasta", package = "MiDAS")
  fasta <- seqinr::read.alignment(fasta_file, format = "fasta")
  fasta <- fasta$seq[[1]]
  expect_equal(paste(hla_alignments[1, ], collapse = ""),
               toupper(fasta)
  ) # check sequence with sequence from fasta

  hla_alignments_trim <- readHlaAlignments(aln_file, trim = TRUE)
  n_trimmed <- ncol(hla_alignments) - ncol(hla_alignments_trim)
  expect_equal(n_trimmed, 24) # check if sequence is trimmed properly

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
               "input file contains no correct HLA alignments"
  )
  unlink(fake_aln_tmp)

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
  expect_error(readHlaAlignments(fake_aln_tmp),
               "start codon is not marked properly in the input file"
  )
  unlink(fake_aln_tmp)
})

test_that("KIR haplotype calls are read properly", {
  file <- system.file("extdata", "KIP_output_example.txt", package = "MiDAS")
  kir_calls <- readKirCalls(file, counts = FALSE)
  kir_calls_test <- read.table(file = file,
                               header = TRUE,
                               sep = "\t",
                               na.strings = c("", "NA"),
                               stringsAsFactors = FALSE
  )
  expect_equal(kir_calls, kir_calls_test)

  kir_calls <- readKirCalls(file)
  kir_counts_test <- kirHaplotypeToCounts(kir_calls_test[, 2, drop = TRUE])
  kir_counts_test <- cbind(ID = kir_calls$ID, kir_counts_test[, -1], stringsAsFactors = FALSE)
  rownames(kir_counts_test) <- NULL
  expect_equal(kir_calls, kir_counts_test)

  expect_error(readKirCalls(file = "foo"), "Path 'foo' does not exist")

  expect_error(readKirCalls(file, hap_dict = "foo"), "Path 'foo' does not exist")

  expect_error(readKirCalls(file, counts = "foo"),
               "counts is not a flag \\(a length one logical vector\\).")

  expect_error(readKirCalls(file, binary = "foo"),
               "binary is not a flag \\(a length one logical vector\\).")

  expect_error(readKirCalls(file, na.strings = 1),
               "na.strings is not a character vector")

  extracol_file <- tempfile()
  kir_calls <- readKirCalls(file, counts = FALSE)
  write.table(kir_calls[, c(1, 2, 2)], file = extracol_file, sep = "\t")
  expect_error(readKirCalls(extracol_file),
               "KIR haplotypes calls table should have 2 columns, not 3")
  unlink(extracol_file)

  badhap_file <- tempfile()
  kir_calls <- readKirCalls(file, counts = FALSE)
  kir_calls[1, 2] <- "foo"
  write.table(kir_calls, file = badhap_file, sep = "\t")
  expect_error(readKirCalls(badhap_file),
               "rows 1 of input file contains unexpected characters")
  unlink(badhap_file)
})
