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
})

test_that("HLA allele resolution is number of sets of digits * 2", {
  expect_equal(getAlleleResolution(c("A*01", "A*01:24", "B*01:25:22",
                                     "C*05:24:55:54")
               ),
               c(2, 4, 6, 8)
  )
})

test_that("Reduced HLA allele have desired resoulution", {
  expect_equal(reduceAlleleResolution(c("A*01", "A*01:24", "B*01:25:22",
                                        "C*05:24:55:54"), 2
               ),
               c("A*01", "A*01", "B*01", "C*05")
  )
})
