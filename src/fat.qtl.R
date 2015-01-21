#####################################################################
# This script runs a QTL analysis of body fat phenotypes from five
#  overlapping subcongenics on MMU2 and the founder congenic.  Each
#  congenic was used to develop an F2 intercross and all the data
#  was merged as a single cross. All mice were first genotyped with
#  a panel of 11 microsattelite markers.  The recombinant mice were
#  later genotyped with a high density panel of 48 SNP in between
#  140 to 180 Mb (roughly) to increase the resoluition of the QTL map
#
# To cite this data please use doi:10.5281/zenodo.12793, and the
#  publication.
#####################################################################

# move to analysis
setwd("../analysis/")

# libraries and sources
library(qtl)
source('../src/qtl.R')
source('../src/qtl2.R')

# options
run.permutations = TRUE

# loads cross data from data/ 
load(file="../data/SBC_hg2d.rda")

# creating a sex and sac covariate matrix
sm=cbind(D2sq$pheno$sex,D2sq$pheno$sac)

# analysis of total fat
#  nc = no covariates ; sex = sex as covariate ; sac = sac as covariate ; ss = sex and sac as covariate
#  ss2 = sex and sac as additive and interactive covariates
# This section was the initial analysis of raw phenotypes without a corretion for strain
out.TF.nc <- scanone(D2sq, pheno.col="tf", method="hk")
out.TF.sex <- scanone(D2sq, pheno.col="tf", method="hk", addcovar=D2sq$pheno$sex)
out.TF.sac <- scanone(D2sq, pheno.col="tf", method="hk", addcovar=D2sq$pheno$sac)
out.TF.ss <- scanone(D2sq, pheno.col="tf", method="hk", addcovar=sm)
out.TF.ss2  <- scanone(D2sq, pheno.col="tf", method="hk", addcovar=sm, intcovar=sm)

# Analysis of total fat using the residuals of a general linear model that accounts for the effects of strain.
#  This step was performed in fat.preproc.R
out.TFe.nc <- scanone(D2sq, pheno.col=find.pheno(D2sq, "tf.e"), method="hk")
out.TFe.sex <- scanone(D2sq, pheno.col=find.pheno(D2sq, "tf.e"), method="hk", addcovar=D2sq$pheno$sex)
out.TFe.sac <- scanone(D2sq, pheno.col=find.pheno(D2sq, "tf.e"), method="hk", addcovar=D2sq$pheno$sac)
out.TFe.ss <- scanone(D2sq, pheno.col=find.pheno(D2sq, "tf.e"), method="hk", addcovar=sm)
out.TFe.ss2  <- scanone(D2sq, pheno.col=find.pheno(D2sq, "tf.e"), method="hk", addcovar=sm, intcovar=sm)


## Analysis of individual fat pads using raw data
out.GFP.nc <- scanone(D2sq, pheno.col="gfp", method="hk")
out.RFP.nc <- scanone(D2sq, pheno.col="rfp", method="hk")
out.MFP.nc <- scanone(D2sq, pheno.col="mfp", method="hk")
out.FFP.nc <- scanone(D2sq, pheno.col="ffp", method="hk")

out.GFP.ss <- scanone(D2sq, pheno.col="gfp", method="hk", addcovar=sm)
out.RFP.ss <- scanone(D2sq, pheno.col="rfp", method="hk", addcovar=sm)
out.MFP.ss <- scanone(D2sq, pheno.col="mfp", method="hk", addcovar=sm)
out.FFP.ss <- scanone(D2sq, pheno.col="ffp", method="hk", addcovar=sm)

out.GFP.ss2 <- scanone(D2sq, pheno.col="gfp", method="hk", addcovar=sm, intcovar=sm)
out.RFP.ss2 <- scanone(D2sq, pheno.col="rfp", method="hk", addcovar=sm, intcovar=sm)
out.MFP.ss2 <- scanone(D2sq, pheno.col="mfp", method="hk", addcovar=sm, intcovar=sm)
out.FFP.ss2 <- scanone(D2sq, pheno.col="ffp", method="hk", addcovar=sm, intcovar=sm)


# Analysis of individual fat pads using the residuals from a general linear model that accounts for the effects of strain.
#  This step was performed in fat.preproc.R
out.GFPe.nc <- scanone(D2sq, pheno.col="gfp.e", method="hk")
out.RFPe.nc <- scanone(D2sq, pheno.col="rfp.e", method="hk")
out.MFPe.nc <- scanone(D2sq, pheno.col="mfp.e", method="hk")
out.FFPe.nc <- scanone(D2sq, pheno.col="ffp.e", method="hk")

out.GFPe.ss <- scanone(D2sq, pheno.col="gfp.e", method="hk", addcovar=sm)
out.RFPe.ss <- scanone(D2sq, pheno.col="rfp.e", method="hk", addcovar=sm)
out.MFPe.ss <- scanone(D2sq, pheno.col="mfp.e", method="hk", addcovar=sm)
out.FFPe.ss <- scanone(D2sq, pheno.col="ffp.e", method="hk", addcovar=sm)

out.GFPe.ss2 <- scanone(D2sq, pheno.col="gfp.e", method="hk", addcovar=sm, intcovar=sm)
out.RFPe.ss2 <- scanone(D2sq, pheno.col="rfp.e", method="hk", addcovar=sm, intcovar=sm)
out.MFPe.ss2 <- scanone(D2sq, pheno.col="mfp.e", method="hk", addcovar=sm, intcovar=sm)
out.FFPe.ss2 <- scanone(D2sq, pheno.col="ffp.e", method="hk", addcovar=sm, intcovar=sm)

# Analisis of multiple traits (Compatible with current versions of R/qtl, but no longer as useful)
#  These functions are contained in the qtl.R and qtl2.R.  At the time of the analysis R/qtl did not
#  have built in functions to analyze multiple QTL, in addition these functions allow for searching
#  for multiple peaks in the lod profile. Ploting functions automaitcally plot the CI of each QTL
#  above the peak.
#sm=cbind(D2sq$pheno$sex,D2sq$pheno$sac)
sw= cbind(D2sq$pheno$sex, D2sq$pheno$WK2)

out.FATtrte.nc <- trait.scan(D2sq, c('tf.e', 'gfp.e', 'rfp.e', 'mfp.e', 'ffp.e'), method='hk')
out.FATtrte.ss <- trait.scan(D2sq, c('tf.e', 'gfp.e', 'rfp.e', 'mfp.e', 'ffp.e'), method='hk', addcov=sm)
out.FATtrte.ss2 <- trait.scan(D2sq, c('tf.e', 'gfp.e', 'rfp.e', 'mfp.e', 'ffp.e'), method='hk', addcov=sm, intcov=sm)

out.FATtrt.nc <- trait.scan(D2sq, c('tf', 'gfp', 'rfp', 'mfp', 'ffp'), method='hk')
out.FATtrt.ss <- trait.scan(D2sq, c('tf', 'gfp', 'rfp', 'mfp', 'ffp'), method='hk', addcov=sm)
out.FATtrt.ss2 <- trait.scan(D2sq, c('tf', 'gfp', 'rfp', 'mfp', 'ffp'), method='hk', addcov=sm, intcov=sm)

mpeaks.scanone(out.TFe.ss, sep=1, cutoff=3)
mlodint.scanone(out.TFe.ss, chr=2, sep=1, drop=1.5, cutoff=3)

mlodint.mscanone(out.FATtrt.ss, chr=2, sep=1, drop=1.5, cutoff=1)
mlodint.mscanone(out.FATtrt.nc, chr=2, sep=1, drop=1.5, cutoff=2)


save(list = ls(pattern = "out"), file = "scanone_output.rda")

# permutations
if(run.permutations) {
source("../src/fat.perm.R", verbose = TRUE)
}

# Analyzing only HG2D cross
#  This data is not saved in the rda files.

HG2D <- subset(D2sq, ind=grep("2D$", D2sq$pheno$strain))
SBC <- subset(D2sq, ind=grep("2D[ABCDE]$", D2sq$pheno$strain))
ss <- cbind(HG2D$pheno$sex, HG2D$pheno$sac)
ssb <- cbind(SBC$pheno$sex, SBC$pheno$sac)

Fate.ss.2D <- trait.scan(HG2D, c('tf.e', 'gfp.e', 'rfp.e', 'mfp.e', 'ffp.e'), method='hk', addcov=ss)
Fat.ss.2D <- trait.scan(HG2D, c('tf', 'gfp', 'rfp', 'mfp', 'ffp'), method='hk', addcov=ss)
Fat.ss.2D.em <- trait.scan(HG2D, c('tf', 'gfp', 'rfp', 'mfp', 'ffp'), method='em', addcov=ss)
Fate.ss.SBC <- trait.scan(SBC, c('tf.e', 'gfp.e', 'rfp.e', 'mfp.e', 'ffp.e'), method='hk', addcov=ssb)


