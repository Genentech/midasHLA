% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{runMiDAS}
\alias{runMiDAS}
\title{Run MiDAS statistical analysis}
\usage{
runMiDAS(
  object,
  experiment,
  inheritance_model = NULL,
  conditional = FALSE,
  omnibus = FALSE,
  omnibus_groups_filter = NULL,
  lower_frequency_cutoff = NULL,
  upper_frequency_cutoff = NULL,
  correction = "bonferroni",
  n_correction = NULL,
  exponentiate = FALSE,
  th = 0.05,
  th_adj = TRUE,
  keep = FALSE,
  rss_th = 1e-07
)
}
\arguments{
\item{object}{An existing fit from a model function such as lm, glm and many
others.}

\item{experiment}{String indicating the experiment associated with
\code{object}'s \code{MiDAS} data to use. Valid values includes:
\code{"hla_alleles"}, \code{"hla_aa"}, \code{"hla_g_groups"},
\code{"hla_supertypes"}, \code{"hla_NK_ligands"}, \code{"kir_genes"},
\code{"kir_haplotypes"}, \code{"hla_kir_interactions"},
\code{"hla_divergence"}, \code{"hla_het"}, \code{"hla_custom"},
\code{"kir_custom"}. See \code{\link{prepareMiDAS}} for more information.}

\item{inheritance_model}{String specifying inheritance model to use.
Available choices are \code{"dominant"}, \code{"recessive"},
\code{"additive"}.}

\item{conditional}{Logical flag indicating if conditional analysis should be
performed.}

\item{omnibus}{Logical flag indicating if omnibus test should be used.}

\item{omnibus_groups_filter}{Character vector specifying omnibus groups to
use.}

\item{lower_frequency_cutoff}{Number giving lower frequency threshold.
Numbers greater than 1 are interpreted as the number of feature occurrences,
numbers between 0 and 1 as fractions.}

\item{upper_frequency_cutoff}{Number giving upper frequency threshold.
Numbers greater than 1 are interpreted as the number of feature occurrences,
numbers between 0 and 1 as fractions.}

\item{correction}{String specifying multiple testing correction method. See
details for further information.}

\item{n_correction}{Integer specifying number of comparisons to consider
during multiple testing correction calculations. For Bonferroni correction
it is possible to specify a number lower than the number of comparisons
being made. This is useful in cases when knowledge about the biology or
redundance of alleles reduces the need for correction. For other methods it
must be at least equal to the number of comparisons being made; only set
this (to non-default) when you know what you are doing!}

\item{exponentiate}{Logical flag indicating whether or not to exponentiate
the coefficient estimates. Internally this is passed to
\code{\link[broom]{tidy}}. This is typical for logistic and multinomial
regressions, but a bad idea if there is no log or logit link. Defaults to
FALSE.}

\item{th}{Number specifying threshold for a variable to be considered
significant.}

\item{th_adj}{Logical flag indicating if adjusted p-value should be used as
threshold criteria, otherwise unadjusted p-value is used.}

\item{keep}{Logical flag indicating if the output should be a list of results
resulting from each selection step. Default is to return only the final
result.}

\item{rss_th}{Number specifying residual sum of squares threshold at which
function should stop adding additional variables. As the residual sum of
squares approaches \code{0} the perfect fit is obtained making further
attempts at variable selection nonsense. This behavior can be controlled
using \code{rss_th}.}
}
\value{
Analysis results, depending on the parameters:
  \describe{
    \item{\code{conditional=FALSE, omnibus=FALSE}}{
      Tibble with first column \code{"term"} holding names of tested 
      variables (eg. alleles). Further columns depends on the used 
      model and are determined by associated \code{tidy} function. Generally 
      they will include \code{"estimate"}, \code{"std.error"}, 
      \code{"statistic"}, \code{"p.value"}, \code{"conf.low"}, 
      \code{"conf.high"}, \code{"p.adjusted"}.
    }
    \item{\code{conditional=TRUE, omnibus=FALSE}}{
      Tibble or a list of tibbles, see \code{keep} argument. The first column 
      \code{"term"} hold names of tested variables. Further 
      columns depends on the used model and are determined by associated 
      \code{tidy} function. Generally they will include \code{"estimate"}, 
      \code{"std.error"}, \code{"statistic"}, \code{"p.value"}, 
      \code{"conf.low"}, \code{"conf.high"}, \code{"p.adjusted"}.
    }
    \item{\code{conditional=FALSE, omnibus=TRUE}}{
      Tibble with first column holding names of tested omnibus groups 
      (eg. amino acid positions) and second names of variables in the group
      (eg. residues). Further columns are: \code{"df"} giving difference in 
      degrees of freedom between base and extended model, \code{"statistic"}
      giving Chisq statistic, \code{"p.value"} and \code{"p.adjusted"}.
    }
    \item{\code{conditional=TRUE, omnibus=TRUE}}{
      Tibble or a list of tibbles, see \code{keep} argument. The first column 
      hold names of tested omnibus groups (eg. amino acid positions), second 
      column hold names of variables in the group (eg. residues). Further 
      columns are: \code{"df"} giving difference in degrees of freedom 
      between base and extended model, \code{"statistic"} giving Chisq 
      statistic, \code{"p.value"} and \code{"p.adjusted"}.
    }
  }
}
\description{
\code{runMiDAS} perform association analysis on MiDAS data using statistical
model of choice. Function is intended for use with \code{\link{prepareMiDAS}}.
See examples section.
}
\details{
By default statistical analysis is performed iteratively on each variable in
selected experiment. This is done by substituting \code{placeholder} in the
\code{object}'s formula with each variable in the experiment.

Setting \code{conditional} argument to \code{TRUE} will cause the statistical
analysis to be performed in a stepwise conditional testing manner, adding the
previous top-associated variable as a covariate to \code{object}'s formula.
The analysis stops when there is no more significant variables, based on
self-defined threshold (\code{th} argument). Either adjusted or unadjusted
p-values can be used as the selection criteria, which is controlled using
\code{th_adj} argument.

Setting \code{omnibus} argument to \code{TRUE} will cause the statistical
analysis to be performed iteratively on groups of variables (like residues at
particular amino acid position) using likelihood ratio test.

Argument \code{inheritance_model} specifies the inheritance model that should
be applyed to experiment's data. Following choices are available: 
\itemize{
  \item{"dominant"}{
    carrier status is sufficient for expression of the phenotype (non-carrier: 
    0, heterozygous & homozygous carrier: 1).
  }
  \item{"recessive"}{
    two copies are required for expression of the phenotype (non-carrier & 
    heterozygous carrier: 0, homozygous carrier: 1).
  }
  \item{"additive"}{
    allele dosage matters, homozygous carriers show stronger phenotype 
    expression or higher risk than heterozygous carriers (non-carrier = 0, 
    heterozygous carrier = 1, homozygous carrier = 2).
  }
  \item{"overdominant"}{
    heterozygous carriers are at higher risk compared to non-carriers or 
    homozygous carriers (non-carrier & homozygous carrier = 0, heterozygous 
    carrier = 1).
  }
}  
 
\code{correction} specifies p-value adjustment method to use, common choice
is Benjamini & Hochberg (1995) (\code{"BH"}). Internally this is passed to
\link[stats]{p.adjust}.
}
\examples{
# create MiDAS object
midas <- prepareMiDAS(hla_calls = MiDAS_tut_HLA,
                      colData = MiDAS_tut_pheno,
                      experiment = c("hla_alleles", "hla_aa")
)

# construct statistical model
object <- lm(disease ~ term, data = midas)

# run analysis
runMiDAS(object, experiment = "hla_alleles", inheritance_model = "dominant")

# omnibus test
# omnibus_groups_filter argument can be used to restrict omnibus test only
# to selected variables groups, here we restrict the analysis to HLA-A
# positions 29 and 43.
runMiDAS(
  object,
  experiment = "hla_aa",
  inheritance_model = "dominant",
  omnibus = TRUE,
  omnibus_groups_filter = c("A_29", "A_43")
)

}
