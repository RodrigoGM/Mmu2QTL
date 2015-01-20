# QTL Analysis for body fat mass

## 1. fat.preproc.R : Data data pre-processing
This script reads the genotype and phenotype and creates the cross object for R/qtl.  Additionally, we adjust phenotypes for the different covariates as described in the Supplementary methods of the paper.  We proceed to estimate the recombination frequencies, calculate genotye probabilities, estimating genotyping errors, and simulate pseudo-markers.

## 2. fat.qtl.R : QTL analysis of fat phenotypes
This script analyses each fat phenotype of this study.  It processes each phenotype independently and collectively with custom functions.  In addition we performed other analyses not present in the manuscript.

## 3. fat.perm.R : QTL permutations
This script runs the permutations for the selected phenotype models.