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
# Install from Bioconductor
BiocManager::install("midasHLA")

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
protein alignment files for all genes shipped with the package. Those are used
for translating HLA alleles to amino acid level. The alignments are downloaded
from the [IPD-IMGT/HLA GitHub repository](https://github.com/ANHIG/IMGTHLA);
the `release` variable in the script selects the IPD-IMGT/HLA release (eg.
`"3650"` for 3.65.0).

`inst/scripts/parse_alignments.R` script pre-parses alignments files 
for package use. Purpose of using pre-parsed files is to speed up 
allele to amino acid sequence translation. Resulting `.Rdata` files 
should be then placed in `inst/extdata` replacing the old alignment 
files. 

Both scripts load the package with `devtools::load_all()` and use relative
paths, so they have to be run from the root of a clone of this repository:

```
Rscript inst/scripts/download_extdata.R   # downloads files to alignments/
Rscript inst/scripts/parse_alignments.R   # writes *_prot.Rdata to the current directory
mv *_prot.Rdata inst/extdata/
```

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

### Package performance

Package performance can be checked using scripts located in `inst/benchmark`. 
Briefly, we used GNU's `time` (`/usr/bin/time -v`) to measure time and memory 
consumption of data transformation and 2 workflows. All tests include loading 
package and reading in input files steps. Look into `inst/benchmark/*` for more 
details. 

## Citing
Migdal M, Ruan DF, Forrest WF, Horowitz A, Hammer C (2021) MiDAS—Meaningful Immunogenetic Data at Scale. PLOS Computational Biology 17(7): e1009131. https://doi.org/10.1371/journal.pcbi.1009131

