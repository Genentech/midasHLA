## 0.0.9045 - 10/05/20
+ fix bug treating NA values as distinct amino acid, that was adding false columns to amino acids counts

## 0.0.9044 - 10/05/20
+ changes amino acids numbering scheme to addhere with EBI convention where numbering starts from -30 position

## 0.0.9043 - 22/04/20
+ fixes hlaToAAVariation failing on mixed resolution hla_calls input
+ adds warning when allele can not be found in the reference alignment

## 0.0.9042 - 17/04/20
+ removes unused pvalue_cutoff argument from runMiDAS function
+ fixes bug with conditional testing of glm models

## 0.0.9041 - 9/04/20
+ fixes problem with placeholder (dummy variable) being constant. This caused coxme function to throw errors upon model definition.

## 0.0.9040 - 16/02/20
+ add functionality for Grantham distance calcualtion between amino acid sequences (`distGrantham` function) and between HLA alleles (`hlaCallsGranthamDistance`).
+ add Grantham distances calculation for Class I HLA alleles to `prepareMiDAS` and `runMiDAS` functionality.

## 0.0.9039 - 18/12/19
+ now using interactions terms is possible through placeholder argument in runMiDAS and presence of dummy var in input data frame

## 0.0.0.9038 - 19/11/2019
+ now variables passed to runMiDAS via variables argument are not taken into account when applying frequency cut-off

## 0.0.0.9037 - 16/11/2019
+ adds new analysis type "none" to runMiDAS function. In this mode only terms passed by variables argument are used.

## 0.0.0.9036 - 8/10/2019
+ refactor function for kabling results now it is much more general and easier to use.

## 0.0.0.9035 - 6/10/2019
+ reading alignment files shipped with package is now faster.
+ alignment files shipped with the packge are now stored in a preparsed format.
+ readHlaAlignments now have divergent behaviour if file is supplied it behaves as previously, when gene is supplied it loads preparsed alignment ommiting some lengthly steps. This change speeds up the process but has no visible consequences.
+ hlaToAAVariation no longer have alnpath argument. The variable amino acid positions search in a custom alignment file can now be done only manualy.

## 0.0.0.9034 - 1/10/2019
+ add function for models comparison using likelihood ratio test (LRTest)
+ add function for performing omnibus test on amino acid positions (aaPosOmnibusTest)

## 0.0.0.9034 - 1/10/2019
+ in summariseAAPosition output alleles numbers are now sorted alphabeticaly

## 0.0.0.9033 - 30/09/2019
+ imporves documentation
+ renames prepareMiDASData to prepareMiDAS
+ renames analyzeMiDASData to runMiDAS

## 0.0.0.9032 - 21/09/2019
+ add function that lists HLA alleles and amino acid residues at a given position together with frequencies.

## 0.0.0.9031 - 21/09/2019
+ add pattern argument to analyzeMiDASData that can be used to subset variables choosen by analysis_type

## 0.0.0.9030 - 15/09/2019
+ add assert to test if there is appropiate tidy function available

## 0.0.0.9029 - 15/09/2019
+ fixes potential bug where NA is accepted as logical values in is.flag asserts

## 0.0.0.9028 - 15/09/2019
+ fixes bug that was occuring when hla_calls passed to hlaToVariable contained NAs every where except one column.

## 0.0.0.9027 - 15/09/2019
+ automated kabling of results from analyzeMiDASData function has been removed, this can be now done by separate call to formatResults function.

## 0.0.0.9026 - 15/09/2019
+ changes dictionaties naming convention. Now alleles dictionaties are named Match_allele_HLA_name.txt
+ add new functionality for converting counts to variables, new dictionaries are named Match_counts_name.txt

## 0.0.0.9025 - 30/08/2019
+ changes the behavior of hlaToVariable and hlaCallsToCounts. hlaToVariable now labels alleles not present in the dictionary with 0, this can be changed with na.value argument. The NAs that were already present in input hla_calls are preserved. hlaCallsToCounts ignores all values that can be converted to numeric - this way alleles that were not present in dictionary and are represented as 0 are not counted.

## 0.0.0.9024 - 24/08/2019
+ adds n_correction argument to analyzeAssociations, analyzeConditionalAssociations and analyzeMiDAS that can be used to pass n argument to p.adjusted calculations.

## 0.0.0.9023 - 14/08/2019
+ fixes handling NAs when counting variables occurences. Previously presence of NA in any allele was ignored and all not observed alleles were treated as so. Now in cases when both alleles are NA this is reflected in counts. However if only one of the alleles is missing behaviour is still the same as previously. The described new behaviour is spread to hlaToVariable and getHlaKirInteractions functions.

## 0.0.0.9022 - 19/08/2019
+ bahvior of variables argument in analyzeMiDASData is now changed, such that it can be used to supply additional variables to analyse that cannot be selected by analysis_type

## 0.0.0.9021 - 30/07/2019
+ changes the behaviour of analyzeMiDASData conditional = TRUE. Now it return list of results from all iterations and kables only best results from each iteration.

## 0.0.0.9020 - 27/07/2019
+ adds human friendly erros in readHlaCalls

## 0.0.0.9019 - 27/07/2019
+ makes it possible to specify na.strings while reading HLA and KIR calls input files

## 0.0.0.9018 - 27/07/2019
+ fixes a bug when models were accessed and updated in parent frame, making it error prone in more complex examples. Now modeels are evaluated in the enviorment they were created in (one defined in model$terms .Enviorment attribute.

## 0.0.0.9017 - 23/07/2019
+ adds kir_genes analysis type to prepareMiDASData and analyzeMiDASData
+ adds hla_kir_interactions analysis type to prepareMiDASData and analyzeMiDASData

## 0.0.0.9016 - 20/07/2019
+ adds checkKirCountsFormat to assert KIR counts format
+ adds getHlaKirInteractions to get HLA - KIR interactions as new variables

## 0.0.0.9016 - 14/07/2019
+ add new function MiDAS that combines prepareMiDASData and analyzeMiDASData.
+ result returned by MiDAS stores its input hla_calls and transformed data in corresponding attributes allowing
  transformed data to be reused in subsequent analyzes.
  
## 0.0.0.9015 - 10/07/2019
+ adds kirHaplotypeToCounts for converting KIR halplotypes to genes counts
+ adds readKirCalls that allows parsing KIR haplotypes calls output by kpi (https://github.com/droe-nmdp/kpi)

## 0.0.0.9014 - 28/06/2019
+ in analyzeMiDASData now it is possible to specify both lower and upper threshold on frequency, before it was only lower.

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

