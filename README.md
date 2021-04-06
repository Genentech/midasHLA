# Meaningful Immunogenetic Data at Scale <img src="vignettes/MiDAS_logo.png" align="right" alt="" width="120" />

Human immunogenetic variation in the form of HLA and KIR 
types has been shown to be strongly associated with a 
multitude of immune-related phenotypes. We present MiDAS, 
an R package enabling statistical association analysis and 
using immunogenetic data transformation functions for HLA 
amino acid fine mapping, analysis of HLA evolutionary 
divergence as well as HLA-KIR interactions. MiDAS closes the 
gap between inference of immunogenetic variation and its 
efficient utilization to make meaningful discoveries.

## Installation

``` r
# Install development version from GitHub
devtools::install_github("Genentech/MiDAS")
```

## Usage

A user tutorial is available here: 
[https://genentech.github.io/midasHLA/articles/MiDAS_tutorial.html](https://genentech.github.io/MiDAS/articles/MiDAS_tutorial.html)

## Developers notes

### External data

The package is shipped together with external data such as 
alignment files or allele frequencies. With time, it is 
possible that those resources will become outdated. Here 
follows a brief description of scripts that can be used for
updating those data. It should be noted that data storage details,
at the external sources, may be changed. In such circumstances, 
the scripts might become obsolete. Nevertheless, they should be
a good starting point to update those sources in next package 
iterations.

`inst/scripts/download_extdata.R` script is used to download HLA 
alignments files. Those are used for translating HLA alleles to 
amino acid level. The alignments are downloaded from 
[EBI's IPD-IMGT/HLA database](www.ebi.ac.uk/ipd/imgt/hla/).
In some cases alignment files contain sequences for multiple genes,
 those will be split into separate files (eg. DRB genes)

`inst/scripts/parse_alignments.R` script pre-parses alignments files 
for package use. Purpose of using pre-parsed files is to speed up 
allele to amino acid sequence translation. Resulting `.Rdata` files 
should be then placed in `inst/extdata` replacing the old alignment 
files. 

`data-raw/allele_frequencies.R` script fetches HLA allele frequencies 
(for genes A, B, C, DQA1, DQB1, DPA1, DPB1, DRB1, DRB3, DRB4, DRB5) from 
[www.allelefrequencies.net](www.allelefrequencies.net)
database and save it in a usable format in `data` directory. This data is
then available under the `allele_frequencies` variable.

`data-raw/kir_frequencies.R` script fetches KIR genes frequencies 
(for genes 3DL3, 2DS2, 2DL2, 2DL3, 2DP1, 2DL1, 3DP1, 2DL4, 3DL1, 3DS1, 
2DL5, 2DS3, 2DS5, 2DS4, 2DS1, 3DL2) from 
[www.allelefrequencies.net](www.allelefrequencies.net)
database and save it in a usable format in `data` directory. This data is
then available under the `kir_frequencies` variable.
