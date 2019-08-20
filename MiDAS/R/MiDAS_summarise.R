#' Lists HLA alleles for amino acid residues at a given position in the alignment
#'
#' @param hla_calls Character vector containing HLA allele numbers.
#' @genename
#'
#' @return Logical vector of length 1 specifying if \code{allele} follows HLA
#' alleles naming conventions and have desired resolution.
#'
#' @examples
#' allele <- c("A*01:01", "A*01:02")
#' checkAlleleFormat(allele)
#'
#' @importFrom assertthat assert_that
#' @importFrom stringi stri_detect_regex
#' @export

file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
hla_calls <- readHlaCalls(file) %>% select(starts_with("DPB1"))
file <- system.file("extdata", "DPB1_prot.txt", package = "MiDAS")
aln <- readHlaAlignments(file, resolution = 4)
aln <- aln[unlist(hla_calls), ]
varaa <- aln[, 8]
df <- lapply(unique(names(varaa)), function(x) data.frame(aa = varaa[x], hla = x, count = sum(names(varaa) == x)))
df <- do.call(rbind, df)
ggplot(data = df, mapping = aes(x = reorder(hla, -count), y = count, fill = aa)) +
  geom_bar(stat = "identity", position = ) +
  labs(title = "DPB1_AA_8", y = "count", x = "") +
  theme_bw()

# It should somehow come from this
aa_variation <- hlaToAAVariation(hla_calls)

# mby its better to wait with this function until the association test is implemented

### typical workflow ###
library("survival")
hla_calls_file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
hla_calls <- readHlaCalls(hla_calls_file)
pheno_file <- system.file("extdata", "pheno_example.txt", package = "MiDAS")
pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)
covar_file <- system.file("extdata", "covar_example.txt", package = "MiDAS")
covar <- read.table(covar_file, header = TRUE, stringsAsFactors = FALSE)
midas_data <- prepareMiDASData(hla_calls = hla_calls,
                               pheno = pheno,
                               covar = covar,
                               analysis_type = "aa_level",
                               inheritance_model = "additive"
)

object <- lm(OS ~ AGE + SEX, data = midas_data)
res <- analyzeMiDASData(object, analysis_type = "aa_level", kable_output = FALSE)


aa_matrix <- hlaToAAVariation(hla_calls)
f <- function(midas_data,
              aa_pos,
              frequency_cutoff) {
  aa_pos <- paste0(aa_pos, "_")
  aa_df <- select(midas_data, starts_with(aa_pos))
}


