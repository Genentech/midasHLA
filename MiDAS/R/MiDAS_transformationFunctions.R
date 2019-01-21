#' Converts HLA allele numbers to amino acid variation
#'
#' Converts HLA allele numbers data frame to a matrix holding information on
#' amino acid level variation.
#'
#' @param hla_calls Data frame containing HLA allele calls, in a format as
#'                  return by `readHlaCalls` function.
#' @param indels Logical indicating whether indels should be considered as
#'               variability.
#' @param unkchar Logical indicating whether unknown characters in the alignment
#'                should be treated as variability.
#' @param alnpath String providing optional path to directory containing
#'                HLA alignment files. Each alignment file have to be named
#'                following EBI database convention GENENAME_prot.txt. If
#'                \code{alnpath} is provided alignment files shipped with the
#'                package are ignored.
#'
#' @return Matrix containing variable amino acid positions. Rownames corresponds
#'         to ID column of input data frame, and colnames to alignment positions
#'         for given genes. If no variation in amino acids alignments is found
#'         function return one column matrix filled with `NA`.
#'
#' @examples
#' hla_calls <- system.file("extdata/HLAHD_output_example.txt", package = "MiDAS")
#' hla_calls <- readHlaCalls(hla_calls)
#' aa_variation <- hlaToAAVariation(hla_calls)
#'
#' @importFrom assertthat assert_that see_if is.dir is.flag
#' @importFrom stringi stri_split_fixed
#' @export
hlaToAAVariation <- function(hla_calls,
                             indels = TRUE,
                             unkchar = FALSE,
                             alnpath = system.file("extdata", package = "MiDAS")){
  assert_that(
    is.data.frame(hla_calls),
    see_if(nrow(hla_calls) >= 1 & ncol(hla_calls) >= 2,
           msg = "input data frame have to have at least 1 rows and 2 columns"
    ),
    see_if(! all(checkAlleleFormat(hla_calls[, 1]), na.rm = TRUE),
           msg = "first column of input data frame should specify samples id"
    ),
    see_if(all(checkAlleleFormat(unlist(hla_calls[, -1])), na.rm = TRUE),
           msg = "values in input data frame doesn't follow HLA numbers specification"
    ),
    is.flag(indels),
    is.flag(unkchar),
    is.dir(alnpath)
  )
  ids <- hla_calls[, 1]
  hla_calls <- hla_calls[, -1]

  # get names of genes and corresponding resolutions
  gene_names <- vapply(X = colnames(hla_calls),
                       FUN = function(x) stri_split_fixed(x, "_")[[1]][1],
                       FUN.VALUE = character(length = 1)
  )
  gene_names_uniq <- unique(gene_names)

  # discard genes for which no alignment files are available
  availbable_genes <- list.files(
    path = system.file("extdata", package = "MiDAS"),
    pattern = "_prot.txt$"
  )
  availbable_genes <- vapply(
    X = stri_split_fixed(availbable_genes, "_prot.txt"),
    `[[`, 1,
    FUN.VALUE = character(length = 1)
  )
  gene_names_uniq <- gene_names_uniq[gene_names_uniq %in% availbable_genes]

  hla_resolution <- vapply(X = gene_names_uniq,
                           FUN = function(x) {
                             x_numbers <- unlist(hla_calls[, gene_names == x])
                             x_res <- getAlleleResolution(na.omit(x_numbers))
                             return(min(x_res))
                           },
                           FUN.VALUE = numeric(length = 1),
                           USE.NAMES = TRUE
  )

  # read alignment matrices and convert to desired resolution
  hla_aln <- lapply(X = gene_names_uniq,
                    FUN = function(x) {
                      aln_file <- file.path(alnpath, paste0(x, "_prot.txt"))
                      aln <- readHlaAlignments(
                        file = aln_file,
                        resolution = hla_resolution[x],
                        unkchar = "*"
                      )
                      return(aln)
                    }
  )

  # get aa variations for each gene
  aa_variation <- list()
  for (i in 1:length(gene_names_uniq)) {
    x_calls <- hla_calls[, gene_names == gene_names_uniq[i]]

    # check if there is possibility for variability
    x_calls_uniq <- na.omit(unique(unlist(x_calls)))
    if (length(x_calls_uniq) <= 1) next()

    # get variable aa positions
    hla_aln[[i]] <- hla_aln[[i]][x_calls_uniq, ]
    var_pos <- getVariableAAPos(hla_aln[[i]],
                                varchar = sprintf("[A-Z%s%s]",
                                                  ifelse(indels, "\\.", ""),
                                                  ifelse(unkchar, "\\*", "")
                                )
    )
    var_aln <- lapply(colnames(x_calls), function(allele) {
      mask <- 1:nrow(hla_aln[[i]]) # This is tmp solution as NAs in character index gives oob error
      names(mask) <- rownames(hla_aln[[i]])
      x <- hla_aln[[i]][mask[x_calls[, allele]], var_pos, drop = FALSE]
      colnames(x) <- paste0(allele, "_", "AA_", var_pos)
      return(x)
    })
    var_aln <- do.call(cbind, var_aln)
    ord <- as.vector(vapply(1:length(var_pos),
                  function(j) {
                    c(j, j + length(var_pos))
                  },
                  FUN.VALUE = numeric(length = 2)
    ))
    var_aln <- var_aln[, ord]

    aa_variation[[length(aa_variation) + 1]] <- var_aln
  }
  if (length(aa_variation) > 1) {
    aa_variation <- do.call(cbind, aa_variation)
    rownames(aa_variation) <- ids
  } else if (length(aa_variation) == 1) {
    aa_variation <- aa_variation[[1]]
    rownames(aa_variation) <- ids
  } else {
    aa_variation <- matrix(nrow = length(ids))
    rownames(aa_variation) <- ids
  }

  return(aa_variation)
}
