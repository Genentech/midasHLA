context("MiDAS class longtests")

test_that("prepareMiDAS", {
  hla_calls <- reduceHlaCalls(MiDAS_tut_HLA, resolution = 4)
  kir_calls <- MiDAS_tut_KIR
  phenotype <- MiDAS_tut_pheno

  args_c <-
    expand.grid(
      lower_frequency_cutoff = c(0, 0.54),
      upper_frequency_cutoff = c(1, 0.63),
      stringsAsFactors = FALSE
    )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])

    midas <- prepareMiDAS(
      hla_calls = hla_calls,
      kir_calls = kir_calls,
      colData = phenotype,
      experiment = c(
        "hla_alleles",
        # "hla_aa",
        "hla_g_groups",
        "hla_supertypes",
        "hla_NK_ligands",
        "kir_genes",
        "hla_kir_interactions",
        "hla_divergence"
      ),
      lower_frequency_cutoff = args$lower_frequency_cutoff,
      upper_frequency_cutoff = args$upper_frequency_cutoff
    )
    validObject(midas)
  }

  expect_error(
    prepareMiDAS(),
    "hla_calls or kir_calls argument has to be specified"
  )


  expect_error(
    prepareMiDAS(hla_calls = cars[1:2, ]),
    "values: 2, 10 in hla_calls doesn't follow HLA numbers specification"
  )

  expect_error(
    prepareMiDAS(hla_calls = hla_calls, colData = cars[1:2, ]),
    "first column in colData must be named 'ID'"
  )

  expect_error(
    prepareMiDAS(
      hla_calls = hla_calls,
      colData = phenotype,
      experiment = 1
    ),
    "experiment is not a character vector"
  )

  expect_error(
    prepareMiDAS(
      hla_calls = hla_calls,
      colData = phenotype,
      experiment = "foo"
    ),
    "experiment should match values \"hla_alleles\", \"hla_aa\", \"hla_g_groups\", \"hla_supertypes\", \"hla_NK_ligands\", \"kir_genes\", \"kir_haplotypes\", \"hla_kir_interactions\", \"hla_divergence\", \"hla_het\", \"hla_custom\", \"kir_custom\"."
  )

  expect_error(
    prepareMiDAS(
      hla_calls = hla_calls,
      colData = phenotype,
      experiment = "hla_alleles",
      placeholder = 1
    ),
    "placeholder is not a string \\(a length one character vector\\)."
  )

  expect_error(
    prepareMiDAS(
      hla_calls = hla_calls,
      colData = phenotype,
      experiment = "hla_alleles",
      placeholder = "outcome"
    ),
    "Placeholder 'outcome' can not be used, it is already used as column name in one of the inputs."
  )
})

test_that("prepareMiDAS_hla_aa", {
  args_c <- expand.grid(
    indels = c(TRUE, FALSE),
    unkchar = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  for (i in 1:nrow(args_c)) {
    args <- as.list(args_c[i, , drop = FALSE])
    args$hla_calls <- MiDAS_tut_HLA
    experiment <- do.call(prepareMiDAS_hla_aa, args)

    counts <-
      hlaToAAVariation(
        hla_calls = args$hla_calls,
        indels = args$indels,
        unkchar = args$unkchar
      )
    counts <-
      aaVariationToCounts(counts)
    counts <- dfToExperimentMat(counts)
    pos <- gsub("_.{1}$", "", rownames(counts))
    pos <- unique(pos)
    omnibus_groups <- lapply(pos, function(p) {
      grep(pattern = paste0("^", p, "_"),
           x = rownames(counts),
           value = TRUE)
    })
    names(omnibus_groups) <- pos
    experiment_test <-
      SummarizedExperiment(assays = counts,
                           metadata = list(
                             inheritance_model_applicable = TRUE,
                             pop_mul = 2,
                             omnibus_groups = omnibus_groups)
                           )
    expect_equal(experiment, experiment_test)
  }
})
