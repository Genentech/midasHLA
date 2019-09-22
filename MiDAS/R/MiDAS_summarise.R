#' Summarize amino acid position
#'
#' Lists HLA alleles and amino acid residues at a given position
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams hlaToAAVariation
#' @param pos String specifying gene and amio acid position, example
#'   \code{"A_9"}.
#' @param aln Matrix containing amino acids sequence alignments as returned by
#'   \link{readHlaAlignments} function. By default function will use alignment
#'   files shipped with the package.
#' @param na.rm Logical flag indicating if \code{NA} values should be considered
#'   for frequency calculations.
#'
#' @return Data frame containing HLA alleles, thier corresponding amino acid
#'   residues and thier frequencies at requested position.
#'
#' @examples
#' file <- system.file("extdata", "HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(file)
#' summariseAAPosition(hla_calls, "A_9")
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr group_by n summarise
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @export
summariseAAPosition <- function(hla_calls,
                                pos,
                                aln = NULL,
                                na.rm = FALSE) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    is.string(pos),
    see_if(grepl("^[A-Z]+[0-9]*_[0-9]+$", pos),
           msg = "amino acid position should be formatted like: A_9."
    ),
    see_if(
      is.matrix(aln) || is.null(aln),
      msg = "aln is not a matrix or NULL."
    ),
    isTRUEorFALSE(na.rm)
  )

  gene <- gsub("[0-9]+$", "", pos)
  pos <- as.integer(gsub(".*_([0-9]+)$", "\\1", pos))

  alleles <- dplyr::select(hla_calls, dplyr::starts_with(gene)) %>%
    unlist()
  alleles_wo_na <- na.omit(alleles)
  assert_that(length(alleles_wo_na) != 0,
              msg = "hla_calls for given gene contains only NA."
  )

  if (na.rm) {
    alleles <- alleles_wo_na
  }

  hla_resolution <- min(getAlleleResolution(alleles_wo_na))
  aln <- readHlaAlignments(
    gene = gsub("_", "", gene),
    resolution = hla_resolution,
    unkchar = "*"
  )

  assert_that(
    see_if(
      pos <= ncol(aln),
      msg = sprintf("amino acid position %i is higher than amino acid sequence length.", pos)
    ),
    see_if(
      all(i <- alleles_wo_na %in% rownames(aln)),
      msg = sprintf(
        fmt = "allele %s could not be found in the nucleotide alignment file.",
        paste(unique(alleles_wo_na[! i]), collapse = ", ")
      )
    )
  )

  aa <- data.frame(
    allele = alleles,
    residue = aln[alleles_wo_na, pos],
    stringsAsFactors = FALSE
  ) %>%
    group_by(.data$allele, .data$residue) %>%
    summarise(count = n(), frequency = .data$count / length(alleles)) %>%
    as.data.frame()

  return(aa)
}
