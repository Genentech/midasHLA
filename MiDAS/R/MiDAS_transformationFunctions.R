#' Converts HLA allele numbers to amino acid variation
#'
#' Converts HLA allele numbers data frame to a data frame holding information on
#' amino acid level varation.
#'
#' @param hla_calls Data frame containing HLA allele calls.
#'
#' @return Matrix containing all variable amino acid positions.
#'
#' @importFrom stringi stri_split_fixed
#' @export
hlaToAAVariation <- function(hla_calls){
  # check if hla_calls is hla_calls
  ids <- hla_calls[, 1]
  hla_calls <- hla_calls[, -1]

  # get names of genes and corresponding resolutions
  gene_names <- vapply(X = colnames(hla_calls),
                       FUN = function(x) stri_split_fixed(x, "_")[[1]][1],
                       FUN.VALUE = character(length = 1)
  )
  hla_resolution <- vapply(X = unique(gene_names),
                           FUN = function(x) {
                             x_numbers <- unlist(hla_calls[, gene_names == x])
                             x_res <- getAlleleResolution(na.omit(x_numbers))
                             return(min(x_res))
                           },
                           FUN.VALUE = numeric(length = 1),
                           USE.NAMES = TRUE
  )

  # read alignment matrices and convert to desired resolution -- if we will switch to SQL db this could be the only chunk changed
  hla_aln <- lapply(X = unique(gene_names),
                    FUN = function(x) {
                      path <- system.file("extdata/",
                                          paste0(x, "_prot.txt"),
                                          package = "MiDAS"
                      )
                      aln <- readHlaAlignments(path)
                      alleles <- rownames(aln)
                      alleles <- reduceAlleleResolution(alleles,
                                                  resolution = hla_resolution[x] # How to format it nicely
                      )
                      unique_idx <- ! duplicated(alleles)
                      aln <- aln[unique_idx, ]
                      rownames(aln) <- alleles[unique_idx]
                      return(aln)
                    }
  )

  # get aa variations for each gene
  aa_variation <- list()
  gene_names_uniq <- unique(gene_names)
  for (i in 1:length(gene_names_uniq)) {
    x_calls <- hla_calls[, gene_names == gene_names_uniq[i]]

    # check if there is possiblity for variability
    x_calls_uniq <- na.omit(unique(unlist(x_calls)))
    if (length(x_calls_uniq) <= 1) next()

    # get variable aa positions
    hla_aln[[i]] <- hla_aln[[i]][x_calls_uniq, ]
    aa_var_pos <- getVariableAAPos(hla_aln[[i]])
    aa_var <- lapply(colnames(x_calls), function(allele) {
      x <- hla_aln[[i]][x_calls[, allele], aa_var_pos, drop = FALSE]
      colnames(x) <- paste0(allele, "_", "AA_", aa_var_pos)
      return(x)
    })
    aa_var <- do.call(cbind, aa_var)
    ord <- as.vector(vapply(1:length(aa_var_pos),
                  function(j) {
                    c(j, j+length(aa_var_pos))
                  },
                  FUN.VALUE = numeric(length = 2)
    ))
    aa_var <- aa_var[, ord]

    aa_variation[[length(aa_variation) + 1]] <- aa_var
  }
  aa_variation <- do.call(cbind, aa_variation)
  rownames(aa_variation) <- ids

  return(aa_variation)
}
