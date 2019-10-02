#!/usr/bin/env R
library(dplyr)
library(broom)

# Proposal of using likelihood ratio test as omnibus test for selecting amino
# acid positions likely associated with phenotype for futher testing.

# likelihood ratio test https://rpubs.com/roes7096/LTR
lrt <- function(mod, mod0) {
  mod_ll <- logLik(mod)
  mod0_ll <- logLik(mod0)
  df <- length(coef(mod)) - length(coef(mod0))
  teststat <- 2 * (as.numeric(mod_ll) - as.numeric(mod0_ll))
  pchisq(teststat, df = df, lower.tail = FALSE)
}

# create data sets
aa_pos <- structure(list(ID = c("PAT1", "PAT2", "PAT3", "PAT4", "PAT5", "PAT6",
                                "PAT7", "PAT8", "PAT9", "PAT10", "PAT11",
                                "PAT12", "PAT13", "PAT14", "PAT15", "PAT16",
                                "PAT17", "PAT18", "PAT19", "PAT20"),
                         A_56_G = c(0L, 2L, 2L, 0L, 2L, 0L, 0L, 0L, 2L, 0L, 2L,
                                    2L, 2L, 2L, 2L, 0L, 2L, 0L, 2L, 2L),
                         A_56_R = c(2L, 0L, 0L, 2L, 0L, 2L, 2L, 2L, 0L, 2L, 0L,
                                    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L),
                         A_56_H = c(0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
                                    0L, 0L, 0L, 0L, 2L, 0L, 2L, 0L, 0L),
                         A_74_H = c(2L, 2L, 0L, 2L, 0L, 0L, 0L, 2L, 0L, 2L, 0L,
                                    2L, 2L, 2L, 1L, 0L, 1L, 0L, 1L, 0L),
                         A_74_D = c(0L, 0L, 2L, 0L,  2L, 2L, 2L, 0L, 2L, 0L, 2L,
                                    0L, 0L, 0L, 1L, 2L, 1L, 2L, 1L, 2L),
                         OS = c(280L, 458L, 415L, 211L, 631L, 183L, 52L, 149L,
                                413L, 227L, 392L, 483L, 501L, 323L, 569L, 554L,
                                413L, 161L, 115L, 499L),
                         AGE = c(46L, 37L, 67L, 75L, 46L, 57L, 43L, 24L, 82L,
                                 46L, 23L, 67L, 35L, 24L, 57L, 65L, 36L, 52L,
                                 74L, 34L)),
                    class = "data.frame", row.names = c(NA, -20L))

mod0 <- lm(OS ~ AGE, data = aa_pos)

mod_A_56 <- lm(OS ~ AGE + A_56_G + A_56_R + A_56_H, data = aa_pos)
mod_A_74 <- lm(OS ~ AGE + A_74_H + A_74_D, data = aa_pos)

res_lrt <- data.frame(pos = c("A_56", "A_74"),
                  pvalue = c(lrt(mod_A_56, mod0), lrt(mod_A_74, mod0))
)
res_lrt

# In this example for including position A_56 p-value is below threshold
# (say 0.05) so we will include it in further analysis. For A_74 its not so we
# will exclude it.

# Further analysis in general looks like this
res <- list()
res[[1]] <- lm(OS ~ AGE + A_56_G, data = aa_pos)
res[[2]] <- lm(OS ~ AGE + A_56_R, data = aa_pos)
res[[3]] <- lm(OS ~ AGE + A_56_H, data = aa_pos)

do.call(rbind, lapply(res, tidy)) %>%
  dplyr::filter(grepl("A_[0-9]+_[A-Z]", term))
