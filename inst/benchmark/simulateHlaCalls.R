#!/usr/bin/env R
# Migdal
# Simulate hla_calls data

devtools::load_all()

options(scipen = 999)

argv <- commandArgs(trailingOnly=TRUE)
out_prefix <- argv[1] # output prefix eg. inst/benchmark/hla_calls
genes <- argv[2:length(argv)] # genes for which to generate calls
size <- c(100, 1000, 10000, 100000) # generated data size
phomo = 0.108 # homozygote probability

genHlaData <- function(alleles, aprob, size=100, phomo=0.108) {
  al_1 <- sample(x = alleles, size = size, replace = TRUE, prob = aprob)
  al_2 <- ifelse(
    test = rbinom(n = size, size = 1, prob = phomo), 
    yes = al_1,
    no = sample(x = alleles, size = size, replace = TRUE, prob = aprob))
  df <- data.frame(paste0("SAM", 1:size), al_1, al_2, stringsAsFactors = FALSE)
  gn <- gsub("\\*.*", "", alleles[1])
  colnames(df) <- c("ID", paste0(gn, "_1"), paste0(gn, "_2"))
  
  return(df)
}

simHlaCalls <- function(allele_list, size, den, phomo) {
  out <- list()
  genes <- names(allele_list)
  for (gene in genes) {
    alleles <- allele_list[[gene]]
    aprob <- sample(x = den$x, size = length(alleles), replace = TRUE, prob = den$y)
    out[[length(out) + 1]] <- genHlaData(alleles, aprob, size, phomo)
  }
  out <- Reduce(f=function(...) left_join(..., by = "ID"), x=out)
  
  return(out)
}

# estimate frequencies density based on published data
dis <- filter(allele_frequencies, grepl("A\\*", var)) %>%
  filter(grepl("Poland pop 3", population))
den <- density(dis$frequency, adjust=3, from=0, to=1)

# Plot allele frequency distribution and fitted kernel
# hist(dis$frequency)
# lines(den, lty="dotted")

# get alleles
allele_list <- lapply(
  X = genes, 
  FUN = function (x) readHlaAlignments(gene = x) %>% rownames())
names(allele_list) <- genes

# simulate HLA data
for (s in size) {
  out <- simHlaCalls(allele_list, s, den, phomo)
  write.table(out, file = paste0(out_prefix, "_", s, ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}

