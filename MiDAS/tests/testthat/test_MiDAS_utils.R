context("HLA allele numbers")

test_that("HLA allele numbers have proper format", {
  expect_equal(checkAlleleFormat(c("A*01", "A*01:24", "B*01:25:22",
                                   "C*05:24:55:54")
                                 ), c(TRUE, TRUE, TRUE, TRUE)
               )
  expect_equal(checkAlleleFormat(c("*01", "A*:22", "C*05:24:55:54")
                                 ), c(FALSE, FALSE, FALSE)
              )
})
