context("IPD-IMGT/HLA alignment releases")

test_that("getAlignmentsRelease", {
  expect_match(getAlignmentsRelease(), "^[0-9]+\\.[0-9]+\\.[0-9]+$")
})

test_that("alignmentsReleaseBranch", {
  expect_equal(alignmentsReleaseBranch("3.65.0"), "3650")
  expect_equal(alignmentsReleaseBranch("3.10.0"), "3100")
  expect_error(
    alignmentsReleaseBranch("3.65"),
    "release should be formatted like: 3.65.0"
  )
})

test_that("prepareHlaAlignment", {
  file <- test_path("alignments", "C_prot.txt")
  cached_aln_obj <- prepareHlaAlignment(file)
  aln <- cached_aln_obj[[1]]
  parsed <- readHlaAlignments(file, unkchar = "*")

  expect_equal(aln[rownames(parsed), ], parsed)
  expect_equal(colnames(aln)[cached_aln_obj[[2]]], "1")
  # lower resolution alleles are inferred
  expect_equal(aln["C*17:03", ], parsed["C*17:03:01:01", ])
  expect_equal(aln["C*17:03:01", ], parsed["C*17:03:01:01", ])

  # alignments of multiple genes can be subset to a single gene
  aln <- prepareHlaAlignment(file, gene = "C")[[1]]
  expect_true(all(startsWith(rownames(aln), "C*")))
  expect_error(
    prepareHlaAlignment(file, gene = "B"),
    "alignment does not contain alleles of gene B"
  )
})

test_that("readHlaAlignments reads other releases on demand", {
  cache_dir <- tempfile()
  old_options <- options(
    midasHLA.cache_dir = cache_dir,
    midasHLA.alignments_release = NULL
  )
  on.exit({
    options(old_options)
    unlink(cache_dir, recursive = TRUE)
  })
  downloads <- 0
  local_mocked_bindings(
    downloadHlaAlignment = function(gene, release, dir) {
      downloads <<- downloads + 1
      file <- file.path(dir, paste0(gene, "_prot.txt"))
      file.copy(test_path("alignments", paste0(gene, "_prot.txt")), file)
      file
    }
  )
  expected <- prepareHlaAlignment(test_path("alignments", "C_prot.txt"))[[1]]

  aln <- readHlaAlignments(gene = "C", release = "3.65.0", unkchar = "*")
  expect_equal(aln, expected)
  expect_equal(downloads, 1)
  expect_true(
    file.exists(file.path(cache_dir, "alignments", "3.65.0", "C_prot.Rdata"))
  )

  # cached release is not downloaded again
  aln <- readHlaAlignments(gene = "C", release = "3.65.0", trim = TRUE,
                           unkchar = "*")
  expect_equal(aln, expected[, which(colnames(expected) == "1"):ncol(expected)])
  expect_equal(downloads, 1)

  # release can be set with an option
  options(midasHLA.alignments_release = "3.65.0")
  expect_equal(readHlaAlignments(gene = "C", unkchar = "*"), expected)
  expect_equal(downloads, 1)

  # shipped release is read from the package
  options(midasHLA.alignments_release = getAlignmentsRelease())
  shipped <- readHlaAlignments(gene = "C")
  options(midasHLA.alignments_release = NULL)
  expect_equal(shipped, readHlaAlignments(gene = "C"))
  expect_equal(downloads, 1)

  expect_error(
    readHlaAlignments(gene = "C", release = "3.65"),
    "release should be formatted like: 3.65.0"
  )
})

test_that("historical IPD-IMGT/HLA releases can be parsed", {
  # downloads alignments, run only when explicitly requested, e.g. on CI
  skip_if_not(
    identical(Sys.getenv("MIDASHLA_TEST_RELEASES"), "true"),
    "set MIDASHLA_TEST_RELEASES=true to test historical releases"
  )
  cache_dir <- tempfile()
  old_options <- options(midasHLA.cache_dir = cache_dir)
  on.exit({
    options(old_options)
    unlink(cache_dir, recursive = TRUE)
  })

  # releases are always downloaded and parsed, also the shipped one, as the
  # test checks compatibility of the parser with formats of the releases
  parseRelease <- function(gene, release) {
    readRDS(releaseAlignmentFile(gene, release))[[1]]
  }
  releases <- c("3.30.0", "3.44.0", "3.55.0", "3.65.0")
  for (gene in c("C", "DRB3")) {
    alns <- lapply(releases, parseRelease, gene = gene)
    names(alns) <- releases
    for (release in releases) {
      aln <- alns[[release]]
      pos <- colnames(aln)
      numbered <- pos[! grepl(".", pos, fixed = TRUE)]
      # numbered positions are consecutive, without position 0
      first <- as.integer(numbered[1])
      expect_equal(
        numbered,
        as.character(setdiff(seq(first, first + length(numbered)), 0)),
        info = paste(gene, release)
      )
    }
    # residues of alleles present in all releases agree at numbered positions
    alleles <- Reduce(intersect, lapply(alns, rownames))
    for (release in releases[-1]) {
      pos <- intersect(colnames(alns[[1]]), colnames(alns[[release]]))
      pos <- pos[! grepl(".", pos, fixed = TRUE)]
      x <- alns[[1]][alleles, pos]
      y <- alns[[release]][alleles, pos]
      known <- ! x %in% c("*", "") & ! y %in% c("*", "")
      expect_gt(mean(x[known] == y[known]), 0.999, label = paste(gene, release))
    }
  }

  # conserved cysteines of HLA class I heavy chain
  aln <- parseRelease("C", "3.30.0")
  expect_equal(
    unname(aln["C*01:02:01:01", c("101", "164", "203", "259")]),
    rep("C", 4)
  )
})
