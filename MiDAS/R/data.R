#' MiDAS tutorial phenotype data
#'
#' Example phenotype data used in MiDAS tutorial
#'
#' @format Data frame with 1000 rows and 4 columns.
#'   \describe{
#'     \item{ID}{Character sample ID}
#'     \item{disease}{Integer}
#'     \item{lab_value}{Numeric}
#'     \item{outcome}{Integer}
#'   }
#'
"MiDAS_tut_pheno"

#' MiDAS tutorial KIR data
#'
#' Example KIRR presence/absence data used in MiDAS tutorial
#'
#' @format Data frame with 1000 rows and 17 columns. First column holds samples
#'   ID's, following columns holds presence/absence indicators for different
#'   KIR genes.
#'   \describe{
#'     \item{ID}{Character sample ID}
#'     \item{KIR3DL3}{Integer}
#'     \item{KIR2DS2}{Integer}
#'     \item{KIR2DL2}{Integer}
#'     \item{KIR2DL3}{Integer}
#'     \item{KIR2DP1}{Integer}
#'     \item{KIR2DL1}{Integer}
#'     \item{KIR3DP1}{Integer}
#'     \item{KIR2DL4}{Integer}
#'     \item{KIR3DL1}{Integer}
#'     \item{KIR3DS1}{Integer}
#'     \item{KIR2DL5}{Integer}
#'     \item{KIR2DS3}{Integer}
#'     \item{KIR2DS5}{Integer}
#'     \item{KIR2DS4}{Integer}
#'     \item{KIR2DS1}{Integer}
#'     \item{KIR3DL2}{Integer}
#'   }
#'
"MiDAS_tut_KIR"

#' MiDAS tutorial HLA data
#'
#' Example HLA calls data used in MiDAS tutorial
#'
#' @format Data frame with 1000 rows and 19 columns. First column holds samples
#'   ID's, following columns holds HLA alleles calls for different genes.
#'   \describe{
#'     \item{ID}{Character sample ID}
#'     \item{A_1}{Character}
#'     \item{A_2}{Character}
#'     \item{B_1}{Character}
#'     \item{B_2}{Character}
#'     \item{C_1}{Character}
#'     \item{C_2}{Character}
#'     \item{DPA1_1}{Character}
#'     \item{DPA1_2}{Character}
#'     \item{DPB1_1}{Character}
#'     \item{DPB1_2}{Character}
#'     \item{DQA1_1}{Character}
#'     \item{DQA1_2}{Character}
#'     \item{DQB1_1}{Character}
#'     \item{DQB1_2}{Character}
#'     \item{DRA_1}{Character}
#'     \item{DRA_2}{Character}
#'     \item{DRB1_1}{Character}
#'     \item{DRB1_2}{Character}
#'   }
#'
"MiDAS_tut_HLA"

#' Alleles frequencies scraped from allelefrequencies.net
#'
#' Accessed on 28.07.20
#'
#' A dataset containing allele frequencies across 5697 alleles
#' For details visit the search results page in the allelefrequencies.net
#' database website.
#'
#' @format A data frame with 2096 rows and 3 variables:
#' \describe{
#'   \item{var}{allele number, character}
#'   \item{population}{reference population name, character}
#'   \item{frequency}{allele frequency in reference population, float}
#' }
#'
#' @source \url{www.allelefrequencies.net}
"allele_frequencies"

#' Grantham distance
#'
#' Integer vector giving Grantham distance values between pairs of amino acid
#' residues.
#'
#' @format Named integer vector of length 400.
#'
"dict_dist_grantham"

#' KIR genes frequencies scraped from allelefrequencies.net
#'
#' Accessed on 28.08.20
#'
#' A dataset containing KIR genes frequencies across 16 genes.
#' For details visit the search results page in the allelefrequencies.net
#' database website.
#'
#' @format A data frame with 3744 rows and 3 variables:
#' \describe{
#'   \item{var}{allele number, character}
#'   \item{population}{reference population name, character}
#'   \item{frequency}{KIR genes carrier frequency in reference population, float}
#' }
#'
#' @source \url{www.allelefrequencies.net}
"kir_frequencies"
