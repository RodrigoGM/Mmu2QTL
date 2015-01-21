## QTL Analysis for body fat mass

#### 1. fat.preproc.R : Data pre-processing

This script reads the genotype and phenotype and creates the cross object for R/qtl.  Additionally, we adjust phenotypes for the different covariates as described in the Supplementary methods of the paper.  We proceed to estimate the recombination frequencies, calculate genotye probabilities, estimating genotyping errors, and simulate pseudo-markers.

#### 2. fat.qtl.R : IM QTL analysis of fat phenotypes

This script analyzes each fat phenotype of this study.  It processes each phenotype independently and collectively with custom functions.  In addition we performed other analyses not present in the manuscript.

#### 3. fat.perm.R : QTL permutations

** Requires `sm` matrix from fat.qtl.R

It performs 1000 permutations for the selected phenotype models.

#### 4. fat.cim.R : CIM QTL analysis fo fat phenotypes

** Requires object `out.TFe.ss`. Can be found in 'scanone_output.rda' or run in the same session as fat.qtl.R

This script analyzes Total Fat Mass with Composite Interval Mapping.  It runs several instances of CIM with various combination of parameters and covariates

#### 5. fat.f2sim.R : F2 intercross simulation

This script simulates an F2 cross segragating two QTL on the same strand spaced 25 cM appart.  The cross consists of three phenotypes with variance, additive, and dominance values equivalent to Total Fat Mass, 1920 individuals, and two marker spacing scenarios. Cross is analyzed in the same script by IM, CIM and replicated CIM

#### 6. fat.rcim.preproc.R : Reading data and pre-processing

This script functions as script 1. fat.preproc.R, however, it performs additional sim.geno to calculate four step sizes for composite interval mapping.

#### 7-11. fat.rcim_s0.R, fat.rcim_s025.R, fat.rcim_s05.R, fat.rcim_s1.R

These scripts run N ```times``` replications of the function ```cim``` and summarize the results into a single object.  The parameter ```window``` of the ```cim()``` function is modified systematically to 1 cM, 0.5 cM, and 0.25 cM.

#### 12. fat.rcim_smry_s0,1,05.R

This script summarizes Total Fat Mass replicates into ```scanone``` ojbects for integration with R/qtl and other scripts here.

#### 13. fat.plots.R : Exploratory and summary figures

This script generates a collection of figures for QC and interpretation of the results.


## function libraries

#### qtl.R

Contains several functions that extend older versions of R/qtl.  It allows for running 'scanone' on multiple traits at the same time, generating an object class 'mscanone'.  Additional functions have been developed for the 'mscanone' class, e.g. multiple peak summaries and plotting

#### qtl2.R

Contains base function for running replicated composite interval mapping.

#### qtl_Ext.R

Depracated.
