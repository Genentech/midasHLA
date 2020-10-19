## ----setup, echo=FALSE, warning=FALSE, include=FALSE--------------------------
knitr::opts_chunk$set(echo = TRUE, fig.align="center")
devtools::load_all()
library(ggplot2)
library(ggpubr)
library(cowplot)

## ----get_data, echo=TRUE, warning=FALSE, include=TRUE-------------------------
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

## ----show_imported_data , echo=FALSE, warning=FALSE, include=FALSE------------
kable(dat_HLA, caption = "HLA data as imported by MiDAS") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----get_frequencies, echo=TRUE, warning=FALSE--------------------------------
freq_HLA <-
  getHlaFrequencies(dat_HLA,
  compare = TRUE
) %>% 
  filter(Freq > 0.05)

freq_KIR <- getKIRFrequencies(MiDAS_tut_KIR)

## ----show_frequencies, echo=FALSE, warning=FALSE------------------------------
kable(freq_HLA, caption = "HLA frequencies, compared to published references") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

kable(freq_KIR, caption = "KIR frequencies") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----plot_frequencies, echo=FALSE, warning=FALSE------------------------------
plot_HLAfreq <- ggdotchart(freq_HLA, x = "allele", y = "Freq",
           add = "segments",
           sorting = "descending",
           add.params = list(color = "lightgray", size = 1.5),
           rotate = TRUE,
           color = "#0073C2FF",
           palette = "jco",
           dot.size = 2,
           ggtheme = theme_pubr(),
           ylab = "Allele Frequency (%)",
           position = position_dodge(0.9)) + 
    theme_pubclean() +  
    scale_y_continuous(limits = c(0, 1), labels = formattable::percent, breaks = seq(0, 1, by = .2))

plot_HLAfreq

## ----prep_HLA, echo=TRUE, warning=FALSE---------------------------------------
HLA <- prepareMiDAS(
  hla_calls = dat_HLA,
  colData = dat_pheno,
  experiment = c("hla_alleles", "hla_NK_ligands")
)

## ----run_HLA_allelic, echo=TRUE, warning=FALSE--------------------------------
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

## ----show_HLA_allelic, echo=FALSE, warning=FALSE------------------------------
kableResults(HLA_results) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----run_HLA_allelic_cond, echo=TRUE, warning=FALSE---------------------------
HLA_results_cond <- runMiDAS(
  object = HLA_model, 
  experiment = "hla_alleles", 
  conditional = TRUE,
  lower_frequency_cutoff = 0.02, 
  upper_frequency_cutoff = 0.98, 
  correction = "bonferroni", 
  exponentiate = TRUE)

## ----show_HLA_allelic_cond, echo=FALSE, warning=FALSE-------------------------
kableResults(HLA_results_cond) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----run_KIR_ligand, echo=TRUE, warning=FALSE---------------------------------
KIRlig_model <- glm(disease ~ term, data = HLA, family = binomial())
KIRlig_results <- runMiDAS(
  KIRlig_model,
  experiment = "hla_NK_ligands",
  lower_frequency_cutoff = 0.02,
  upper_frequency_cutoff = 0.98,
  correction = "bonferroni",
  exponentiate = TRUE
  )

## ----show_KIR_ligand_results, echo=FALSE, warning=FALSE-----------------------
kableResults(KIRlig_results)

## ----run_KIR_ligand_cond, echo=TRUE, warning=FALSE----------------------------
KIRlig_cond_results <- runMiDAS(
  KIRlig_model,
  experiment = "hla_NK_ligands",
  lower_frequency_cutoff = 0.02,
  upper_frequency_cutoff = 0.98,
  correction = "bonferroni",
  exponentiate = TRUE,
  conditional = TRUE
  )

## ----show_KIR_ligand_cond_results, echo=FALSE, warning=FALSE------------------
kableResults(KIRlig_cond_results) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----run_HLA_AA, echo=TRUE, warning=FALSE-------------------------------------
HLA_AA <- prepareMiDAS(
  hla_calls = dat_HLA,
  colData = dat_pheno,
  experiment = "hla_aa"
  )

HLA_AA_model <- glm(disease ~ term, data = HLA_AA, family = binomial())
HLA_AA_results <- runMiDAS(
  HLA_AA_model,
  experiment = "hla_aa",
  lower_frequency_cutoff = 0.02,
  upper_frequency_cutoff = 0.98,
  correction = "bonferroni",
  exponentiate = TRUE
  )

## ----show_HLA_AA, echo=FALSE, warning=FALSE-----------------------------------
kableResults(HLA_AA_results) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----get_KIR, echo=TRUE, warning=FALSE----------------------------------------
dat_KIR <- readKirCalls(file = system.file("extdata", "MiDAS_tut_KIR.txt", package = "MiDAS"))

## ----show_KIR, echo=FALSE, warning=FALSE--------------------------------------
kable(dat_KIR, caption = "KIR data as imported by MiDAS") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----KIR_frequency, echo=TRUE, warning=FALSE----------------------------------
kir_freq <- getKIRFrequencies(dat_KIR) %>% 
  filter(Freq > 0.05)

## ----show_KIR_frequency, echo=FALSE, warning=FALSE----------------------------
kable(kir_freq , caption = "KIR data as imported by MiDAS") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----run_KIR_gene, echo=TRUE, warning=FALSE-----------------------------------
HLAKIR <- prepareMiDAS(
  hla_calls = dat_HLA,
  kir_calls = dat_KIR,
  colData = dat_pheno,
  experiment = c("kir_genes", "hla_kir_interactions")
  )
HLAKIR_model <- glm(disease ~ term, data = HLAKIR, family = binomial())
KIR_results <- runMiDAS(
  HLAKIR_model,
  experiment = "kir_genes",
  lower_frequency_cutoff = 0.02,
  upper_frequency_cutoff = 0.98,
  exponentiate = TRUE
  )

## ----show_KIR_results, echo=FALSE, warning=FALSE------------------------------
kableResults(KIR_results) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----HLA-KIR_interactions_frequency-------------------------------------------
kir_freq <- getFrequencies(HLAKIR, experiment = "hla_kir_interactions")

## ----show_HLA-KIR_interactions_frequency, echo=FALSE, warning=FALSE-----------
kable(kir_freq, caption = "HLA-KIR interactions frequencies") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----run_HLA-KIR_interactions, echo=TRUE, warning=FALSE-----------------------
HLA_KIR_results <- runMiDAS(
  HLAKIR_model,
  experiment = "hla_kir_interactions",
  lower_frequency_cutoff = 0.02,
  upper_frequency_cutoff = 0.98,
  exponentiate = TRUE
  )

## ----show_HLA-KIR_results, echo=FALSE, warning=FALSE--------------------------
kableResults(HLA_KIR_results) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

## ----run_HLAhet, echo=TRUE, warning=FALSE-------------------------------------
### Add heterozygosity as analysis_type, once possible. Ideally, one could analyze both at the same time: analysis_type = c("hla_het","hla_divergence").

HLAdiv <- prepareMiDAS(hla_calls = dat_HLA, colData = dat_pheno, experiment = c("hla_divergence"))
HLAdiv_model <- glm(outcome ~ term, data=HLAdiv, family=binomial())
HLAdiv_results <-
  runMiDAS(HLAdiv_model, experiment = "hla_divergence", exponentiate = TRUE)

## ----show_HLAhet, echo=FALSE, warning=FALSE-----------------------------------
kableResults(HLAdiv_results) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
  scroll_box(width = "100%", height = "200px")

