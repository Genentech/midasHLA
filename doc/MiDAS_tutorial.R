## ----setup, echo=FALSE, warning=FALSE, include=FALSE--------------------------
knitr::opts_chunk$set(echo = TRUE, fig.align="center")
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(kableExtra)
library(MiDAS)
library(tidyr)

## ----logo, echo=FALSE, warning=FALSE------------------------------------------
htmltools::img(
  src = knitr::image_uri("MiDAS_logo.png"),
  alt = 'logo',
  style = 'position:absolute; top:0; right:0; padding:10px; width:250px'
  )

## ----get_data, echo=TRUE, warning=FALSE---------------------------------------
dat_pheno <-
  read.table(
  file = system.file("extdata", "MiDAS_tut_pheno.txt", package = "MiDAS"),
  header = TRUE
  )
dat_HLA <-
  readHlaCalls(
  file = system.file("extdata", "MiDAS_tut_HLA.txt", package = "MiDAS"),
  resolution = 4
  )

## ----show_imported_data, echo=FALSE, warning=FALSE----------------------------
dat_HLA %>%
  head(20) %>%
  kable(caption = "HLA data as imported by MiDAS") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "300px")

## ----get_frequencies, echo=TRUE, warning=FALSE--------------------------------
freq_HLA <- getHlaFrequencies(hla_calls = dat_HLA, compare = TRUE) %>%
  filter(Freq > 0.01)

## ----show_frequencies, echo=FALSE, warning=FALSE------------------------------
freq_HLA %>%
  head() %>%
  kable(caption = "HLA frequencies, compared to published references") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

## ----plot_frequencies, echo=TRUE, warning=FALSE, fig.height=7, fig.width=7----
freq_HLA_long <- gather(
  data = freq_HLA,
  key = "population",
  value = "freq",
  "Freq",
  "USA NMDP European Caucasian",
  "USA NMDP Chinese",
  "USA NMDP African American pop 2",
  factor_key = TRUE
) %>% 
  filter(grepl("^A", allele))

plot_HLAfreq <-
  ggplot(data = freq_HLA_long, aes(x = allele, y = freq, fill = population)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(0.7),
    width = 0.7,
    colour = "black"
  ) +
  coord_flip() +
  scale_y_continuous(labels = formattable::percent)

plot_HLAfreq

## ----prep_HLA, echo=TRUE, warning=FALSE---------------------------------------
HLA <- prepareMiDAS(
  hla_calls = dat_HLA,
  colData = dat_pheno,
  experiment = "hla_alleles"
)

## ----HWE_filtering, echo=TRUE, warning=FALSE----------------------------------
HLA <- HWETest(
  object = HLA,
  experiment = "hla_alleles",
  HWE_cutoff = 0.05 / 447,
  as.MiDAS = TRUE
)

## ----run_HLA_allelic, echo=TRUE, warning=FALSE, message=FALSE-----------------
HLA_model <- glm(disease ~ term, data = HLA, family = binomial())
HLA_results <- runMiDAS(
  object = HLA_model, 
  experiment = "hla_alleles", 
  inheritance_model = "dominant",
  lower_frequency_cutoff = 0.02, 
  upper_frequency_cutoff = 0.98, 
  correction = "bonferroni", 
  exponentiate = TRUE
)

kableResults(HLA_results)

## ----run_HLA_allelic_cond, echo=TRUE, warning=FALSE, message=FALSE------------
HLA_results_cond <- runMiDAS(
  object = HLA_model, 
  experiment = "hla_alleles", 
  inheritance_model = "dominant", 
  conditional = TRUE,
  lower_frequency_cutoff = 0.02, 
  upper_frequency_cutoff = 0.98, 
  correction = "bonferroni", 
  exponentiate = TRUE
)

kableResults(HLA_results_cond, scroll_box_height = "200px")

## ----prepare_HLA_AA, echo=TRUE, warning=FALSE---------------------------------
HLA_AA <- prepareMiDAS(
  hla_calls = dat_HLA,
  colData = dat_pheno,
  experiment = "hla_aa"
)

## ----HLA_AA_to_df, echo=TRUE, warning=FALSE-----------------------------------
dat_HLA_AA <- HLA_AA[["hla_aa"]] %>% 
  assay() %>% 
  t() %>% 
  as.data.frame() %>% 
  select(starts_with("B_97_")) %>% 
  head()

## ----show_HLA_AA, echo=FALSE, warning=FALSE-----------------------------------
kable(dat_HLA_AA, caption = "HLA amino acid data as inferred by MiDAS") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

## ----run_HLA_AA_omnibus, echo=TRUE, warning=FALSE, message=FALSE--------------
HLA_AA_model <- glm(disease ~ term, data = HLA_AA, family = binomial())
HLA_AA_omnibus_results <- runMiDAS(
  HLA_AA_model,
  experiment = "hla_aa",
  conditional = FALSE,
  omnibus = TRUE,
  lower_frequency_cutoff = 0.02,
  upper_frequency_cutoff = 0.98,
  correction = "bonferroni"
)

kableResults(HLA_AA_omnibus_results)

## ----run_HLA_AA_DQB1_9, echo=TRUE, warning=FALSE------------------------------
HLA_AA_DQB1_9_results <- runMiDAS(
  HLA_AA_model,
  experiment = "hla_aa",
  omnibus_groups_filter = "DQB1_9",
  lower_frequency_cutoff = 0.02,
  upper_frequency_cutoff = 0.98,
  correction = "bonferroni",
  exponentiate = TRUE
)

kableResults(HLA_AA_DQB1_9_results, scroll_box_height = "250px")

## ----alleles_for_aa, echo=TRUE, warning=FALSE---------------------------------
HLA_AA_DQB1_9_alleles <- getAllelesForAA(HLA_AA,"DQB1_9")

## ----get_HLA_AA_DQB1_9_alleles, echo=FALSE, warning=FALSE---------------------
kable(HLA_AA_DQB1_9_alleles, escape = FALSE) %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

## ----prep_KIR_ligand, echo=TRUE, warning=FALSE--------------------------------
NKlig <- prepareMiDAS(
  hla_calls = dat_HLA,
  colData = dat_pheno,
  experiment = c("hla_alleles", "hla_NK_ligands")
)

## ----run_KIR_ligand, echo=TRUE, warning=FALSE, message=FALSE------------------
NKlig_model <- glm(disease ~ term, data = NKlig, family = binomial())
NKlig_results <- runMiDAS(
  NKlig_model,
  experiment = "hla_NK_ligands",
  inheritance_model = "dominant",
  correction = "bonferroni",
  exponentiate = TRUE
)

kableResults(NKlig_results)

## ----get_KIR, echo=TRUE, warning=FALSE----------------------------------------
dat_KIR <- readKirCalls(
  file = system.file("extdata", "MiDAS_tut_KIR.txt", package = "MiDAS")
)

## ----show_KIR, echo=FALSE, warning=FALSE--------------------------------------
kable(dat_KIR, caption = "KIR data as imported by MiDAS") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))  %>%
  scroll_box(height = "400px")

## ----KIR_frequency, echo=TRUE, warning=FALSE----------------------------------
kir_freq <- getKIRFrequencies(dat_KIR)

## ----show_KIR_frequency, echo=FALSE, warning=FALSE----------------------------
kable(kir_freq , caption = "KIR data as imported by MiDAS") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(height = "400px")

## ----run_KIR_gene, echo=TRUE, warning=FALSE, message=FALSE--------------------
HLAKIR <- prepareMiDAS(
  hla_calls = dat_HLA,
  kir_calls = dat_KIR,
  colData = dat_pheno,
  experiment = c("hla_NK_ligands","kir_genes", "hla_kir_interactions")
)
KIR_model <- glm(disease ~ term, data = HLAKIR, family = binomial())
KIR_results <- runMiDAS(
  KIR_model,
  experiment = "kir_genes",
  lower_frequency_cutoff = 0.02,
  upper_frequency_cutoff = 0.98,
  exponentiate = TRUE
)

kableResults(KIR_results)

## ----HLA-KIR_interactions_frequency, echo=TRUE, warning=FALSE-----------------
hlakir_freq <- getFrequencies(HLAKIR, experiment = "hla_kir_interactions")

## ----show_HLA-KIR_interactions_frequency, echo=FALSE, warning=FALSE-----------
kable(hlakir_freq, caption = "HLA-KIR interaction frequencies") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(height = "400px")

## ----run_HLA-KIR_interactions, echo=TRUE, warning=FALSE, message=FALSE--------
HLAKIR_model <- glm(disease ~ term, data = HLAKIR, family = binomial())
HLA_KIR_results <- runMiDAS(
  HLAKIR_model,
  experiment = "hla_kir_interactions",
  lower_frequency_cutoff = 0.02,
  upper_frequency_cutoff = 0.98,
  exponentiate = TRUE
)

kableResults(HLA_KIR_results)

## ----run_HLAhet, echo=TRUE, warning=FALSE, message=FALSE----------------------

HLA_het <- prepareMiDAS(
  hla_calls = dat_HLA, 
  colData = dat_pheno, 
  experiment = c("hla_het","hla_divergence")
)

HLA_het_model <- glm(outcome ~ term, data=HLA_het, family=binomial())

HLA_het_results <- runMiDAS(HLA_het_model, 
  experiment = "hla_het", 
  exponentiate = TRUE
)

kableResults(HLA_het_results)

HLA_div_results <- runMiDAS(HLA_het_model, 
  experiment = "hla_divergence", 
  exponentiate = TRUE
)

kableResults(HLA_div_results, scroll_box_height = "250px")

