## 0.0.0.9008 - 11/05/2019
+ adds script for conversion of HLA calls to VCF format

## 0.0.0.9008 - 11/05/2019
+ introduces countsToHlaCalls that allows converting HLA counts table back to
  HLA calls data frame under additive inheritance model 

## 0.0.0.9007 - 8/05/2019
+ forwardConditionalSelection has been renamed to analyzeConditionalAssociations,
  now it returns tibble containing results from conditional testing as the new
  covariates are beeing added. For each variable results from the test upon which
  the variable is used for the first time are showed.
+ Warnings about uninitialized variables comming from dplyr functions has been 
  solved.

## 0.0.0.9006 - 6/04/2019
+ analyzing associations and stepwise selection now has been rewriten such that
  they operate on model objects returned by functions such as lm. The amount of
  inputs is minimized to only variables that are to be analyzed.
+ hlaAssocModels function has been removed as it became obsolete.
+ hlaToAAVariation now returns data frame by default, matrix still can be
  returned for now if specified.
+ Functions for calculation alleles and amino acids frequencies were added.
+ Function for converting amino acids data frame have been added.


## 0.0.0.9005 - 18/03/2019
+ now hla calls are merged with phenotypic and covariate data before using 
  statistical functions, using prepareHlaData. This workflow is going to be
  further strenghtend in future releases.
+ introduce inheritance models specifying how hla calls should be converted to
  hla counts.
+ Response and covariate variables in statistical functions now have to be typed 
  by hand.

## 0.0.0.9004 - 18/03/2019
+ introduce forwardConditionalSelection for stepwise conditional testing, adding 
  the previous top-associated allele as covariate, until there’s no more significant 
  alleles 

## 0.0.0.9003 - 05/02/2019
+ G groups alleles are now accepted as properly formatted
+ G groups alleles can be now reduced, with warning

## 0.0.0.9002 - 27/02/2019
+ introduce analyzeHlaAssociations for performing statistical analysys of HLA
  alleles associations

## 0.0.0.9001 - 29/01/2019
+ added functions for converting alleles into additional variables 

## 0.0.0.9000 - 29/01/2019
### Begining of the CHANGELOG 

