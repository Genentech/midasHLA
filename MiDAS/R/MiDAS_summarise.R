#' Summarize amino acid position
#'
#' List HLA alleles and amino acid residues at a given position.
#'
#' @inheritParams checkHlaCallsFormat
#' @inheritParams hlaToAAVariation
#' @param aa_pos String specifying gene and amino acid position, example
#'   \code{"A_9"}.
#' @param aln Matrix containing amino acid sequence alignments as returned by
#'   \code{\link{readHlaAlignments}} function. By default function will use
#'   alignment files shipped with the package.
#' @param na.rm Logical flag indicating if \code{NA} values should be considered
#'   for frequency calculations.
#'
#' @return Data frame containing HLA alleles, their corresponding amino acid
#'   residues and frequencies at requested position.
#'
#' @examples
#' summariseAAPosition(MiDAS_tut_HLA, "A_9")
#'
#' @importFrom assertthat assert_that
#' @importFrom dplyr group_by n select starts_with summarise
#' @importFrom formattable percent
#' @importFrom magrittr %>%
#' @importFrom stats na.omit
#' @importFrom rlang .data
#' @export
summariseAAPosition <- function(hla_calls,
                                aa_pos,
                                aln = NULL,
                                na.rm = FALSE) {
  assert_that(
    checkHlaCallsFormat(hla_calls),
    is.string(aa_pos),
    see_if(grepl("^[A-Z]+[0-9]*_[0-9]+$", aa_pos),
           msg = "amino acid position should be formatted like: A_9."
    ),
    see_if(
      is.matrix(aln) || is.null(aln),
      msg = "aln is not a matrix or NULL."
    ),
    isTRUEorFALSE(na.rm)
  )

  gene <- gsub("[0-9]+$", "", aa_pos)
  aa_pos <- gsub(".*_([0-9]+)$", "\\1", aa_pos)

  alleles <- select(hla_calls, starts_with(gene)) %>%
    unlist()
  alleles_wo_na <- alleles[! is.na(alleles)]
  assert_that(length(alleles_wo_na) != 0,
              msg = "hla_calls for given gene contains only NA."
  )

  if (na.rm) {
    alleles <- alleles_wo_na
  }

  hla_resolution <- min(getAlleleResolution(alleles_wo_na))
  alleles_wo_na <- reduceAlleleResolution(alleles_wo_na, hla_resolution)
  aln <- readHlaAlignments(
    gene = gsub("_", "", gene),
    resolution = hla_resolution,
    unkchar = "*"
  )

  assert_that(
    see_if(
      aa_pos %in% colnames(aln),
      msg = sprintf("amino acid position %s was not found in amino acid sequence.", aa_pos)
    ),
    see_if(
      all(i <- alleles_wo_na %in% rownames(aln)),
      msg = sprintf(
        fmt = "allele %s could not be found in the alignment file.",
        paste(unique(alleles_wo_na[! i]), collapse = ", ")
      )
    )
  )

  aa <- data.frame(
    allele = gsub("[A-Z]+[0-9]*", "", alleles_wo_na),
    residue = aln[alleles_wo_na, aa_pos],
    stringsAsFactors = FALSE
  ) %>%
    group_by(.data$residue) %>%
    summarise(
      allele = paste(sort(unique(.data$allele)), collapse = ", "),
      count = n(),
      frequency = percent(.data$count / length(alleles))
    ) %>%
    dplyr::select(
      !! sprintf("HLA-%s (%s)", gsub("_", "", gene), aa_pos) := .data$residue,
      !! sprintf("HLA-%s alleles", gsub("_", "", gene)) := .data$allele,
      .data$count,
      .data$frequency
    ) %>%
    as.data.frame()

  return(aa)
}
