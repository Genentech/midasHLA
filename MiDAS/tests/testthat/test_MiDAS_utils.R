context("HLA allele numbers")

test_that("HLA allele numbers have proper format", {
  expect_equal(checkAlleleFormat(c("A*01", "A*01:24", "B*01:25:22",
                                   "C*050:24:55:54")
               ),
               c(TRUE, TRUE, TRUE, TRUE)
  )
  expect_equal(checkAlleleFormat(c("01", "01:24", "01:25:22",
                                   "050:24:55:54")
               ),
               c(FALSE, FALSE, FALSE, FALSE)
  )
  expect_equal(checkAlleleFormat(c("*01", "A*:22", "C*05:24:55:54:89")
               ),
               c(FALSE, FALSE, FALSE)
  )
  expect_error(checkAlleleFormat(1), "allele is not a character vector")
})

test_that("HLA allele resolution is number of sets of digits * 2", {
  expect_equal(getAlleleResolution(c("A*01", "A*01:24", "B*01:25:22",
                                     "C*05:24:55:54")
               ),
               c(2, 4, 6, 8)
  )
  expect_error(getAlleleResolution("word"),
               "allele have to be a valid HLA allele number"
  )
})

test_that("Reduced HLA allele have desired resoulution", {
  expect_equal(reduceAlleleResolution(c("A*01", "A*01:24", "B*01:25:22",
                                        "C*05:24:55:54"), 2
               ),
               c("A*01", "A*01", "B*01", "C*05")
  )
  expect_error(reduceAlleleResolution("C*05:24:55:54", resolution = "four"),
               "resolution is not a count \\(a single positive integer\\)"
  )
  expect_error(reduceAlleleResolution("word", resolution = 4),
               "input resolution can't be lower than requested resolution"
  )
})

context("HLA allele alignments")

test_that("Variable amino acids positions are detected properly", {
  hlaa_calls <- c("A*01:01", "A*01:02")
  hlaa_res <- 4
  hlaa_aln <- readHlaAlignments(system.file("extdata",
                                            "A_prot.txt",
                                            package = "MiDAS")
  )
  four_dig_numbers <- reduceAlleleResolution(rownames(hlaa_aln), resolution = 4)
  hlaa_aln <- hlaa_aln[!duplicated(four_dig_numbers), ]
  rownames(hlaa_aln) <- four_dig_numbers[!duplicated(four_dig_numbers)]
  hlaa_aln <- hlaa_aln[hlaa_calls, ]

  expect_equal(getVariableAAPos(hlaa_aln), c(9, 17))

  expect_error(getVariableAAPos(hlaa_calls), "alignment is not a matrix")

  hlaa_aln_empty_col <- hlaa_aln
  hlaa_aln_empty_col[, 2] <- rep(".", nrow(hlaa_aln_empty_col))
  expect_error(getVariableAAPos(hlaa_aln_empty_col),
               "alignment contains columns with out any amino acid letter"
  )
})
