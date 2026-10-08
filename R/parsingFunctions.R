#' Read HLA allele calls
#'
#' \code{readHlaCalls} read HLA allele calls from file
#'
#' Input file has to be a tsv formatted table with a header. First column should
#' contain sample IDs, further columns hold HLA allele numbers. See
#' \code{system.file("extdata", "MiDAS_tut_HLA.txt", package = "midasHLA")} file
#' for an example.
#'
#' \code{resolution} parameter can be used to reduce HLA allele numbers. If
#' reduction is not needed \code{resolution} can be set to 8. \code{resolution}
#' parameter can take the following values: 2, 4, 6, 8. For more details
#' about HLA allele numbers resolution see
#' \url{http://hla.alleles.org/nomenclature/naming.html}.
#'
#' @inheritParams reduceAlleleResolution
#' @inheritParams utils::read.table
#' @param file Path to input file.
#'
#' @return HLA calls data frame. First column hold sample IDs, further columns 
#'   hold HLA allele numbers.
#'
#' @examples
#' file <- system.file("extdata", "MiDAS_tut_HLA.txt", package = "midasHLA")
#' hla_calls <- readHlaCalls(file)
#'
#' @importFrom assertthat assert_that is.readable see_if
#' @importFrom stats na.omit
#' @importFrom stringi stri_split_fixed
#' @importFrom utils read.table
#' @export
readHlaCalls <- function(file,
                         resolution = 4,
                         na.strings = c("Not typed", "-", "NA")) {
  assert_that(is.readable(file),
              is.count(resolution),
              is.character(na.strings)
  )
  hla_calls <- read.table(file,
                          header = TRUE,
                          sep = "\t",
                          stringsAsFactors = FALSE,
                          na.strings = na.strings
  )
  assert_that(checkHlaCallsFormat(hla_calls))

  # set colnames based on allele numbers
  gene_names <- vapply(X = 2:ncol(hla_calls),
                       FUN = function(i) {
                         names <- stri_split_fixed(hla_calls[, i], "*")
                         names <- vapply(X = names,
                                FUN = function(x) x[1],
                                FUN.VALUE = character(length = 1)
                         )
                         names <- unique(na.omit(names))
                         assert_that(
                           see_if(length(names) <= 1,
                                  msg = "Gene names in columns are not identical"
                           ),
                           see_if(length(names) != 0,
                                  msg = "One of the columns contains only NA"
                           )
                         )
                         return(names)
                       },
                       FUN.VALUE = character(length = 1)
  )
  ord <- order(gene_names)
  gene_names <- gene_names[ord]
  gene_names <- toupper(gene_names)
  gene_names_id <- unlist(lapply(table(gene_names), seq_len))
  gene_names <- paste(gene_names, gene_names_id, sep = "_")
  hla_calls <- hla_calls[, c(1, ord + 1)]
  colnames(hla_calls) <- c("ID", gene_names)

  hla_calls <- reduceHlaCalls(hla_calls, resolution = resolution)

  return(hla_calls)
}

#' Read HLA allele alignments
#'
#' \code{readHlaAlignments} read HLA allele alignments from file.
#'
#' HLA allele alignment file should follow EBI database format, for details
#' see
#' \url{ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/README.md}.
#'
#' All protein alignment files from the EBI database are shipped with the package.
#' They can be easily accessed using \code{gene} parameter. If \code{gene} is
#' set to \code{NULL}, \code{file} parameter is used instead and alignment is
#' read from the provided file. In EBI database alignments for DRB1, DRB3, DRB4
#' and DRB5 genes are provided as a single file, here they are separated.
#' Additionally, were possible sequences for alleles not present in the alignments have been 
#' inferred based on higher resolution alleles. To this end we have reduced 
#' alleles to 6 and 4 digit resolution and took consensus sequence to represent
#' missing alleles. Positions for which there was no full agreement were marked
#' as unknown.
#' 
#' 
#'
#' @inheritParams readHlaCalls
#' @param gene Character vector of length one specifying the name of a gene for
#'   which alignment is required. See details for further explanations.
#' @param trim Logical indicating if alignment should be trimmed to start codon
#'   of the mature protein.
#' @param unkchar Character to be used to represent positions with unknown
#'   sequence.
#' @param release String giving IPD-IMGT/HLA release (eg. \code{"3.44.0"}) of
#'   the alignment to use when reading alignment for \code{gene}. By default the
#'   release of alignments shipped with the package is used, see
#'   \code{\link{getAlignmentsRelease}}. Alignments of other releases are
#'   downloaded from the IMGTHLA GitHub repository
#'   (\url{https://github.com/ANHIG/IMGTHLA}), parsed and cached in the
#'   directory given by the \code{midasHLA.cache_dir} option, by default
#'   \code{tools::R_user_dir("midasHLA", "cache")}. Setting the
#'   \code{midasHLA.alignments_release} option changes the release used by all
#'   functions using HLA alignments, eg. \code{prepareMiDAS}.
#'
#' @return Matrix containing HLA allele alignments.
#'
#'   Rownames correspond to allele numbers and columns to positions in the
#'   alignment. Positions are numbered according to the IPD-IMGT/HLA
#'   nomenclature: residues of the reference allele are numbered from the first
#'   codon of the mature protein (position 1), positions of the signal peptide
#'   are negative and position 0 is omitted. Columns where the reference allele
#'   has no residue hold insertions present in other alleles, they are named
#'   after the preceding position followed by the insertion index, eg.
#'   \code{"4.1"}, \code{"4.2"}. Column names should be treated as text, not
#'   numbers. Sequences following the termination codon are marked as empty
#'   character (\code{""}). Unknown sequences are marked with a character of
#'   choice, by default \code{""}. Stop codons are represented by a hash (X).
#'   Insertion and deletions are marked with period (.). See
#'   \code{vignette("MiDAS_alignments", package = "midasHLA")} for details.
#'
#' @examples
#' hla_alignments <- readHlaAlignments(gene = "A")
#'
#' @importFrom assertthat assert_that is.count is.readable is.string see_if
#' @importFrom stringi stri_flatten stri_split_regex stri_sub
#' @importFrom stringi stri_subset_fixed stri_subset_regex stri_read_lines
#' @importFrom stringi stri_detect_regex
#' @export
readHlaAlignments <- function(file,
                              gene = NULL,
                              trim = FALSE,
                              unkchar = "",
                              release = getOption("midasHLA.alignments_release")) {
  assert_that(
    isTRUEorFALSE(trim),
    is.string(unkchar),
    see_if(
      is.null(release) || (is.string(release) && grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", release)),
      msg = "release should be formatted like: 3.65.0"
    )
  )

  if (is.null(gene)) {
    assert_that(is.readable(file))
    aln_raw <- stri_read_lines(file)
    aln <- stri_split_regex(aln_raw, "\\s+")

    # extract lines containing alignments and omit empty alignment lines
    nonempty_lines <- vapply(aln, length, integer(length = 1)) >= 2
    aln <- aln[nonempty_lines]
    allele_numbers <- vapply(aln, `[`, character(length = 1), 2)
    allele_lines <- checkAlleleFormat(allele_numbers)

    assert_that(
      see_if(any(allele_lines),
             msg = "could not find alleles numbers in the alignment file"
      )
    )
    allele_numbers <- allele_numbers[allele_lines]
    aln <- aln[allele_lines]
    assert_that(
      see_if(all(stri_detect_regex(unlist(aln), "^[A-Z0-9:*.-]*$")),
                               msg = "alignments lines contain non standard characters"
      )
    )

    tmp_aln_env <- new.env(size = 5000)
    for (i in seq_along(allele_numbers)) {
      assign(
        x = allele_numbers[i],
        value = append(
          x = get0(
            x = allele_numbers[i],
            envir = tmp_aln_env,
            ifnotfound = character(length = 0)
          ),
          values = aln[[i]][-c(1, 2)] # discard empty element and allele number
        ),
        envir = tmp_aln_env
      )
    }
    aln_list <- as.list(tmp_aln_env)[unique(allele_numbers)] # convert to list and sort

    # alignment spans the longest allele, as some alleles extend beyond the
    # reference sequence
    aln_list <- lapply(aln_list, stri_flatten)
    seq_along_aln <- seq_len(max(nchar(unlist(aln_list))))
    ref_seq <- stri_sub(aln_list[[1]],
                        seq_along_aln,
                        seq_along_aln
    )
    aln <- do.call(rbind,
                   lapply(aln_list,
                          function(a) {
                            a <- stri_sub(a,
                                          seq_along_aln,
                                          seq_along_aln
                            )
                            i <- a == "-"
                            a[i] <- ref_seq[i]

                            return(a)
                          }
                   )
    )

    # find AA positions numbers
    aln_raw <- aln_raw[nonempty_lines]
    # line with positions numbers, older releases also have a title line
    # containing word "Protein"
    raw_first_codon_idx <- nchar(stri_subset_regex(aln_raw, "^\\s*Prot\\s")[1])
    raw_alignment_line <- stri_sub(aln_raw[allele_lines][1],
                                   1,
                                   raw_first_codon_idx
    )
    raw_alignment_seq <- stri_split_regex(raw_alignment_line, "\\s+")
    raw_alignment_seq <- unlist(raw_alignment_seq)[-c(1, 2)]
    first_codon_idx <- nchar(stri_flatten(raw_alignment_seq))
    assert_that(
      see_if(is.count(first_codon_idx),
             msg = "start codon is not marked properly in the input file"
      )
    )

    # positions are numbered according to the reference allele, starting from
    # the first codon of the mature protein and omitting 0. Columns where the
    # reference allele has no residue are insertions present in other alleles,
    # they are named after the preceding position, eg. 4.1, 4.2, ...
    ref_residues <- ! ref_seq %in% c(".", "")
    assert_that(
      see_if(isTRUE(ref_residues[first_codon_idx]),
             msg = "start codon is not marked properly in the input file"
      )
    )
    ref_count <- cumsum(ref_residues)
    pos <- ref_count - ref_count[first_codon_idx] + 1
    pos <- ifelse(pos <= 0, pos - 1, pos)
    ins <- seq_along(ref_count) - match(ref_count, ref_count) + 1 -
      (ref_count > 0)
    aln_colnames <- ifelse(ref_residues, pos, paste0(pos, ".", ins))
    colnames(aln) <- aln_colnames

    # discard aa '5 to start codon of mature protein
    if (trim) {
      aln <- aln[, first_codon_idx:ncol(aln)]
    }
  } else {
    assert_that(is.string(gene))
    gene <- toupper(gene)
    file <- paste0(system.file("extdata", package = "midasHLA"),
                   "/",
                   gene,
                   "_prot.Rdata"
    )
    available_genes <- list.files(
      path = system.file("extdata", package = "midasHLA"),
      pattern = "_prot.Rdata$",
      full.names = TRUE
    )
    assert_that(
      see_if(file %in% available_genes,
             msg = sprintf("alignment for %s is not available", gene)
      )
    )
    if (! is.null(release) && release != getAlignmentsRelease()) {
      file <- releaseAlignmentFile(gene, release)
    }

    cached_aln_obj <- readRDS(file) # list(readHlaAlignments(file, trim = FALSE, unkchar = "*"), first_codon_idx)
    aln <- cached_aln_obj[[1]]

    # discard aa '5 to start codon of mature protein
    if (trim) {
      first_codon_idx <- cached_aln_obj[[2]]
      aln <- aln[, first_codon_idx:ncol(aln)]
    }
  }

  # substitute unkchar
  aln[aln == "*"] <- unkchar

  return(aln)
}

#' Get IPD-IMGT/HLA release of shipped HLA alignments
#'
#' \code{getAlignmentsRelease} returns the IPD-IMGT/HLA release of the HLA
#' protein alignments shipped with the package.
#'
#' Alignments of other IPD-IMGT/HLA releases can be used by setting the
#' \code{midasHLA.alignments_release} option, see
#' \code{\link{readHlaAlignments}}.
#'
#' @return String giving IPD-IMGT/HLA release, eg. \code{"3.65.0"}.
#'
#' @examples
#' getAlignmentsRelease()
#'
#' @export
getAlignmentsRelease <- function() {
  file <- system.file("extdata", "alignments_release.txt", package = "midasHLA")
  release <- readLines(file, n = 1, warn = FALSE)

  return(trimws(release))
}

#' Get IPD-IMGT/HLA repository branch of a release
#'
#' The IMGTHLA GitHub repository keeps a branch for each IPD-IMGT/HLA release,
#' named after the release number without dots, eg. \code{"3650"} for release
#' \code{"3.65.0"}.
#'
#' @param release String giving IPD-IMGT/HLA release, eg. \code{"3.65.0"}.
#'
#' @return String giving name of the release branch.
#'
#' @importFrom assertthat assert_that is.string see_if
#'
alignmentsReleaseBranch <- function(release) {
  assert_that(
    see_if(
      is.string(release) && grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", release),
      msg = "release should be formatted like: 3.65.0"
    )
  )

  return(gsub(".", "", release, fixed = TRUE))
}

#' Download HLA protein alignment of an IPD-IMGT/HLA release
#'
#' Downloads HLA protein alignment file of a given IPD-IMGT/HLA release from
#' the IMGTHLA GitHub repository (\url{https://github.com/ANHIG/IMGTHLA}).
#' Releases prior to 3.5x provide alignments of DRB genes in a single file
#' (\code{DRB_prot.txt}), which is downloaded instead of a missing gene's
#' file.
#'
#' @param gene String giving name of HLA gene.
#' @inheritParams alignmentsReleaseBranch
#' @param dir String giving path to the directory where the file is saved.
#'
#' @return String giving path to the downloaded file.
#'
#' @importFrom assertthat assert_that see_if
#' @importFrom utils download.file
#'
downloadHlaAlignment <- function(gene, release, dir) {
  branch <- alignmentsReleaseBranch(release)
  url <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/%s/alignments/%s_prot.txt"
  old_options <- options(timeout = max(600, getOption("timeout")))
  on.exit(options(old_options))

  download <- function(name) {
    file <- file.path(dir, paste0(name, "_prot.txt"))
    ok <- tryCatch(
      download.file(sprintf(url, branch, name), destfile = file, quiet = TRUE) == 0,
      error = function(e) FALSE,
      warning = function(w) FALSE
    )
    if (ok) file else NULL
  }
  file <- download(gene)
  if (is.null(file) && grepl("^DRB[0-9]$", gene)) {
    file <- download("DRB")
  }
  assert_that(
    see_if(
      ! is.null(file),
      msg = sprintf(
        "alignment for %s could not be downloaded from IPD-IMGT/HLA release %s",
        gene,
        release
      )
    )
  )

  # check that the file comes from the requested release
  header <- readLines(file, n = 10, warn = FALSE)
  release_pattern <- paste0("IPD-IMGT/HLA.* ", gsub(".", "\\.", release, fixed = TRUE), "$")
  assert_that(
    see_if(
      any(grepl(release_pattern, header)),
      msg = sprintf("downloaded alignment is not from IPD-IMGT/HLA release %s", release)
    )
  )

  return(file)
}

#' Prepare HLA alignment for package use
#'
#' \code{prepareHlaAlignment} reads HLA protein alignment file and infers
#' sequences of lower resolution alleles not present in the alignment. Alleles
#' are reduced to 6 and 4 digit resolution and consensus sequence is used to
#' represent missing alleles. Unknown residues of partially sequenced alleles
#' are ignored; positions where the known residues disagree, or no residue is
#' known, are marked as unknown (\code{"*"}).
#'
#' @inheritParams readHlaAlignments
#' @param gene String giving name of HLA gene. If specified, only alleles of
#'   this gene are kept, which is used for alignment files containing multiple
#'   genes.
#'
#' @return List with two elements: matrix containing HLA allele alignments, as
#'   returned by \code{\link{readHlaAlignments}} with \code{unkchar = "*"}, and
#'   index of the column holding position 1. This is the format of the
#'   alignments shipped with the package.
#'
#' @importFrom assertthat assert_that see_if
#'
prepareHlaAlignment <- function(file, gene = NULL) {
  alignment <- readHlaAlignments(file, trim = FALSE, unkchar = "*")

  if (! is.null(gene)) {
    alignment <- alignment[startsWith(rownames(alignment), paste0(gene, "*")), , drop = FALSE]
    assert_that(
      see_if(
        nrow(alignment) > 0,
        msg = sprintf("alignment does not contain alleles of gene %s", gene)
      )
    )
  }

  # infer missing lower resolution alleles
  for (res in c(6, 4)) {
    allele_numbers <- reduceAlleleResolution(rownames(alignment), resolution = res)
    missing_alleles <- unique(allele_numbers[! allele_numbers %in% rownames(alignment)])
    missing_aln <- lapply(missing_alleles, function(allele) {
      i <- allele_numbers == allele
      apply(alignment[i, , drop = FALSE], 2, function(col) {
        # unknown residues of partially sequenced alleles are ignored
        known <- unique(col[col != "*"])
        if (length(known) == 1) known else "*"
      })
    })
    if (length(missing_aln) == 0) next
    missing_aln <- do.call(rbind, missing_aln)
    rownames(missing_aln) <- missing_alleles
    alignment <- rbind(alignment, missing_aln)
  }

  first_codon_idx <- which(colnames(alignment) == "1")

  return(list(alignment, first_codon_idx))
}

#' Get path to HLA alignment of an IPD-IMGT/HLA release
#'
#' Returns path to the prepared HLA protein alignment of a given IPD-IMGT/HLA
#' release. Alignments are downloaded and prepared once and stored in the cache
#' directory, given by the \code{midasHLA.cache_dir} option, by default
#' \code{tools::R_user_dir("midasHLA", "cache")}.
#'
#' @inheritParams downloadHlaAlignment
#'
#' @return String giving path to the prepared alignment file.
#'
releaseAlignmentFile <- function(gene, release) {
  cache_dir <- getOption(
    "midasHLA.cache_dir",
    default = tools::R_user_dir("midasHLA", "cache")
  )
  release_dir <- file.path(cache_dir, "alignments", release)
  file <- file.path(release_dir, paste0(gene, "_prot.Rdata"))
  if (! file.exists(file)) {
    message(sprintf(
      "Downloading and parsing HLA-%s alignment from IPD-IMGT/HLA release %s",
      gene,
      release
    ))
    download_dir <- tempfile()
    dir.create(download_dir)
    on.exit(unlink(download_dir, recursive = TRUE))
    aln_file <- downloadHlaAlignment(gene, release, download_dir)
    cached_aln_obj <- prepareHlaAlignment(aln_file, gene = gene)
    dir.create(release_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(cached_aln_obj, file = file)
  }

  return(file)
}

#' Read KIR calls
#'
#' \code{readKirCalls} read KIR calls from file.
#'
#' Input file has to be a tsv formatted table. First column should be named
#' "ID" and contain samples IDs, further columns should hold KIR genes presence
#' / absence indicators. See
#' \code{system.file("extdata", "MiDAS_tut_KIR", package = "midasHLA")} for an
#' example.
#'
#' @inheritParams utils::read.table
#' @param file Path to input file.
#'
#' @return Data frame containing KIR gene's counts. First column hold samples 
#'   IDs, further columns hold KIR genes presence / absence indicators.
#'
#' @examples
#' file <- system.file("extdata", "MiDAS_tut_KIR.txt", package = "midasHLA")
#' readKirCalls(file)
#'
#' @importFrom assertthat assert_that is.readable see_if
#' @importFrom dplyr left_join select
#' @importFrom stats na.omit setNames
#' @export
readKirCalls <- function(file,
                         na.strings = c("", "NA", "uninterpretable")) {
  assert_that(
    is.readable(file),
    is.character(na.strings)
  )

  kir_calls <- read.table(
    file = file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    na.strings = na.strings
  )
  checkKirCallsFormat(kir_calls)

  return(kir_calls)
}
