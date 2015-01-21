#################################################################
# This script runs the permutations of the IM QTL analysis of
#  body fat phenotypes from five overlapping subcongenics on MMU2
#  and the founder congenic.  
#
# To cite this data please use doi:10.5281/zenodo.12793, and the
#  publication.
#################################################################
setwd("../analysis/")
library(qtl)

## source from within fat.qtl.R

# permutation for residuals
out.TFe.nc.p <- scanone(D2sq, pheno.col="tf.e", method="hk", n.perm=1000, verbose=TRUE)
out.TFe.ss.p <- scanone(D2sq, pheno.col="tf.e", method="hk", addcovar=sm, n.perm=1000, verbose=TRUE)
out.TFe.ss2.p <- scanone(D2sq, pheno.col="tf.e", method="hk", addcovar=sm, intcovar=sm, n.perm=1000, verbose=TRUE)

summary(out.TFe.nc.p)
summary(out.TFe.ss.p)
summary(out.TFe.ss2.p)

layout(mat = matrix(c(1:3,0), nrow = 2, byrow = TRUE))
plot(out.TFe.nc.p, main = "TFe - No Covariates\n1000x")
plot(out.TFe.ss.p, main = "TFe - Sex & Sac as additive\n1000x")
plot(out.TFe.ss2.p, main = "TFe - Sex and Sac as additive and interactive\n1000x")

out.GFPe.p<- scanone(D2sq, pheno.col="gfp.e", method="hk", addcovar=sm, n.perm=1000, verbose=TRUE)
summary(out.GFPe.p)
#LOD thresholds (1000 permutations)
#       lod
#  5%  5.12
#  10% 4.71

out.RFPe.p<- scanone(D2sq, pheno.col="rfp.e", method="hk", addcovar=sm, n.perm=1000, verbose=TRUE)
summary(out.RFPe.p)
#LOD thresholds (1000 permutations)
#       lod
#  5%  5.81
#  10% 5.06

out.MFPe.p<- scanone(D2sq, pheno.col="mfp.e", method="hk", addcovar=sm, n.perm=1000, verbose=TRUE)
summary(out.MFPe.p)
#LOD thresholds (1000 permutations)
#       lod
#  5%  5.00
#  10% 4.51

out.FFPe.p<- scanone(D2sq, pheno.col="ffp.e", method="hk", addcovar=sm, n.perm=1000, verbose=TRUE)
summary(out.FFPe.p)
#LOD thresholds (1000 permutations)
#       lod
#  5%  5.15
#  10% 4.69

apply(D2sq$pheno[,c("gfp.e", "rfp.e", "mfp.e", "ffp.e", "tf.e")],2, function(m) glm(m ~ sex + sac, family=gaussian, data=D2sq$pheno))
apply(D2sq$pheno[,c("gfp.e", "rfp.e", "mfp.e", "ffp.e", "tf.e")],2, function(m) glm(m ~ sex + sac + sex:sac, family=gaussian, data=D2sq$pheno))


# permutations with sex & sac

out.TF.p <- scanone(D2sq, pheno.col="tf", method="hk", addcovar=sm, intcovar=sm, n.perm=1000, verbose=TRUE)
summary(out.TF.p)
#LOD thresholds (1000 permutations)
#       lod
#   5%  5.10
#  10%  4.58

out.GFP.p<- scanone(D2sq, pheno.col="gfp", method="hk", addcovar=sm, intcovar=sm, n.perm=1000, verbose=TRUE)
summary(out.GFP.p)
#LOD thresholds (1000 permutations)
#       lod
#  5%  5.12
#  10% 4.71

out.RFP.p<- scanone(D2sq, pheno.col="rfp", method="hk", addcovar=sm, intcovar=sm, n.perm=1000, verbose=TRUE)
summary(out.RFP.p)
#LOD thresholds (1000 permutations)
#       lod
#  5%  5.81
#  10% 5.06

out.MFP.p<- scanone(D2sq, pheno.col="mfp", method="hk", addcovar=sm, intcovar=sm, n.perm=1000, verbose=TRUE)
summary(out.MFP.p)
#LOD thresholds (1000 permutations)
#       lod
#  5%  5.00
#  10% 4.51

out.FFP.p<- scanone(D2sq, pheno.col="ffp", method="hk", addcovar=sm, intcovar=sm, n.perm=1000, verbose=TRUE)
summary(out.FFP.p)
#LOD thresholds (1000 permutations)
#       lod
#  5%  5.15
#  10% 4.69

ls(pattern='\\.p$')
save(list=ls(pattern='\\.p$'), file="permutations.rda")
