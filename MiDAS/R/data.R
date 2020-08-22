#' MiDAS tutorial HLA data
#'
#' Example HLA calls data used in MiDAS tutorial
#'
#' @format Data frame with 1000 rows and 19 columns. First column holds samples
#'   ID's, following columns holds HLA alleles calls for different genes.
#'
"MiDAS_tut_HLA"

#' MiDAS tutorial KIR data
#'
#' Example KIRR presence/absence data used in MiDAS tutorial
#'
#' @format Data frame with 1000 rows and 17 columns. First column holds samples
#'   ID's, following columns holds presence/absence indicators for different
#'   KIR genes.
#'
"MiDAS_tut_KIR"

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

#' Alleles frequencies scraped from allelefrequencies.net in wide format
#'
#' Accessed on 28.07.20
#'
#' A dataset containing allele frequencies across 5697 alleles
#' For quick yet more details I recommend the search results page in the
#' allelefrequencies.net database website. Data are gathered in wide format.
#' The general reference populations are defined based on different sources. For exampe
#' in case of HLA-A for Hispanic population frequencies from USA NMDP Hispanic South
#' or Central American are used, however in case of HLA-DQA1 it is Chile Santiago pop 2.
#' See \code{allele_frequencies_reference_populations}.
#'
#' @format A data frame with 2096 rows and 8 variables:
#' \describe{
#'   \item{var}{allele number, character}
#'   \item{population}{reference population name, character}
#'   \item{frequency}{allele frequency in reference population, float}
#' }
#'
#' @source \url{www.allelefrequencies.net}
"allele_frequencies"

#' Reference populations used in allele_frequencies
#'
#' A dataset containing names of reference populations used for
#' extracting alleles frequency from allelefrequencies.net. For
#' some of the HLA genes different populations were used.
#'
#' @format A data frame with 2096 rows and 8 variables:
#' \describe{
#'   \item{population}{general population name, character}
#'   \item{A}{names of the source populations used as reference for HLA gene A, character}
#'   \item{B}{names of the source populations used as reference for HLA gene B, character}
#'   \item{C}{names of the source populations used as reference for HLA gene C, character}
#'   \item{DQA1}{names of the source populations used as reference for HLA gene DQA1, character}
#'   \item{DQB1}{names of the source populations used as reference for HLA gene DQB1, character}
#'   \item{DPA1}{names of the source populations used as reference for HLA gene DPA1, character}
#'   \item{DPB1}{names of the source populations used as reference for HLA gene DPB1, character}
#'   \item{DRB1}{names of the source populations used as reference for HLA gene DRB1, character}
#'   \item{DRB3}{names of the source populations used as reference for HLA gene DRB3, character}
#'   \item{DRB4}{names of the source populations used as reference for HLA gene DRB4, character}
#'   \item{DRB5}{names of the source populations used as reference for HLA gene DRB5, character}
#' }
#'
#' @source \url{www.allelefrequencies.net}
"allele_frequencies_reference_populations"

#' Grantham distance
"dict_dist_grantham"
