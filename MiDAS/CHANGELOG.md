## 0.0.0.9013 - 28/06/2019
+ prepareMiDASData now can accept multiple analysis types at once which allows to cerate data frames with all possible variables in one go
+ variables in data frame produced by prepareMiDASData are now labeled with corresponding analysis type
+ analyzeMiDASData now under variables = NULL selects appropaite variables based on labels associated with midas_data

## 0.0.0.9012 - 25/06/2019
+ fixes analyzeMiDASData function crashing on float variables. Now frequency calculations are done only for properly labeled variables on midas_data.

## 0.0.0.9011 - 21/06/2019
+ fixes allele_group preparation scheme, where before after conversion to groups 
  there were no conversion to counts. Additionally now this scheme have been broken 
  into three: allele_g_group, allele_supertypes, allele_group (Bw4/6, C1/2).

## 0.0.0.9010 - 18/06/2019
+ introduces prepareMiDASData which transforms hla_calls according to predefined schemas, like
  amino acid variation, expression level etc.

## 0.0.0.9009 - 9/06/2019
+ adds script for conversion of HLA calls to VCF format
+ adds script for conversion of VCF to HLA calls

## 0.0.0.9008 - 11/05/2019
+ introduces countsToHlaCalls that allows converting HLA counts table back to
  HLA calls data frame under additive inheritance model 
+ introduces analyzeMiDASData which is a higher level abstraction of analyzeConditionalAssociations 
  and analyzeAssociations. It also pretty format results to kabled html or latex.
+ introduces formatResults that allows to format results to kabled html ir latex.

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

