#####################################################################
# This script runs a QTL analysis of body fat phenotypes from five
#  overlapping subcongenics on MMU2 and the founder congenic.  Each
#  congenic was used to develop an F2 intercross and all the data
#  was merged as a single cross. All mice were first genotyped with
#  a panel of 11 microsattelite markers.  The recombinant mice were
#  later genotyped with a high density panel of 48 SNP in between
#  140 to 180 Mb (roughly) to increase the resoluition of the QTL map
# To cite this data please use doi:10.5281/zenodo.12793, and the
#  publication.
#####################################################################
#################################################################
#   The 07-12-10 *.RDA data set contains the data from the 5  ##
#   subcongenics and the data from the HG2D cross             ##
#  The code here will be used to fine map the critical region ##
#  of Fatq2 using a replicated CIM approacn                   ##
#################################################################

# move to analysis
setwd("../analysis/")

# libraries and sources
library(qtl)
library(coda)
library(boa)
library(hdrcde)
source('../src/qtl.R')
source('../src/qtl2.R')
 
# loads cross data from data/ 
load(file="../data/SBC_hg2d.rda")

# creating a sex and sac covariate matrix
sm=cbind(D2sq$pheno$sex,D2sq$pheno$sac)

# analysis of total fat 
# Composite interval mapping uing ss residuals. These were developed in the fat.preproc.R.
#  TF was corrected for strain, sex, and sacrifice weight.
#  nc = no covariates ; sex = sex as covariate ; sac = sac as covariate ; ss = sex and sac as covariate
#  ss2 = sex and sac as additive and interactive covariates
cim.TFss.w2.imp<- cim(D2sq, pheno.col=find.pheno(D2sq, "tf.ss"), n.marcovar=3, window=2, method="hk", imp.method="imp")
cim.TFss.w2.arg<- cim(D2sq, pheno.col=find.pheno(D2sq, "tf.ss"), n.marcovar=3, window=2, method="hk", imp.method="argmax")

# Replicated composite interval mapping.  We run cim 'n' times with the same parameters.  The output is the object containing the
#  median LOD profile of 'n' replicates, and the rcim.results is a matrix containing the LOD profile for each relicate
# 400x
rcim.TFss.w2.imp.400x<- r.cim(D2sq, pheno.col=find.pheno(D2sq, "tf.ss"), times=400, n.marcovar=3, window=2, method="hk", imp.method="imp", save=TRUE)
ALL.TFss.w2.imp.400x<-rcim.results
rcim.TFss.w2.arg.400x<- r.cim(D2sq, pheno.col=find.pheno(D2sq, "tf.ss"), times=400, n.marcovar=3, window=2, method="hk", imp.method="argmax", save=TRUE)
ALL.TFss.w2.arg.400x<-rcim.results

# 50x
rcim.TFss.w2.imp.50x<- r.cim(D2sq, pheno.col=find.pheno(D2sq, "tf.ss"), times=50, n.marcovar=3, window=2, method="hk", imp.method="imp")
rcim.TFss.w2.arg.50x<- r.cim(D2sq, pheno.col=find.pheno(D2sq, "tf.ss"), times=50, n.marcovar=3, window=2, method="hk", imp.method="argmax")

# confidence interval estimates
bayesint(cim.TFss.w2.arg, chr=2)
bayesint(out.TFe.ss, chr=2)
bayesint(rcim.TFss.w2.imp.400x, chr=2)

TFss.rcim<- r.cim(D2sq, pheno.col=find.pheno(D2sq, "tf.ss"), times=100, method="hk", n.marcovar=3, window=1, save=TRUE)
markers<-rownames(rcim.results)

save(list = ls(pattern = "cim"), file = "cim_output.rda")

# evaluation of cim replicates

markers<- rownames(ALL.TFss.w2.arg.400x)

g.pvals.arg<- sapply(markers, function(m) {mc<- mcmc(data=as.matrix(t(ALL.TFss.w2.arg.400x[m,3:402])))
                                          gd<- geweke.diag(mc, frac=.1, frac2=.5)
                                          p.val<- 1-(pnorm(abs(gd$z))-0.5)*2
                                          return(p.val)
                                          })

g.pvals.imp<- sapply(markers, function(m) {mc<- mcmc(data=as.matrix(t(ALL.TFss.w2.imp.400x[m,3:402])))
                                          gd<- geweke.diag(mc, frac=.1, frac2=.5)
                                          p.val<- 1-(pnorm(abs(gd$z))-0.5)*2
                                          return(p.val)
                                          })
                                          
pvals<- matrix(rbind(g.pvals.arg, g.pvals.imp))
                             

# histograms
pdf(file="../figures/Marker histograms 400x.pdf")
sapply(markers, function(m) {layout(matrix(c(1,2,1,2), 2, 2, byrow=TRUE), c(2,2), respect=TRUE)
                            par(mar=c(4,4,4,0))
                            hist(as.numeric(ALL.TFss.w2.arg.400x[m,3:402]), xlab=m, main="Histogram TF LOD m=ARG", breaks=250)
                            par(mar=c(4,2,4,0))
                            hist(as.numeric(ALL.TFss.w2.imp.400x[m,3:402]), xlab=m, main="Histogram TF LOD m=IMP", breaks=250)
                            })
dev.off()

## DENSITY PLOTS
pdf(file="../figures/Marker HDR density 400x.pdf")
sapply(markers, function(m) {layout(matrix(c(1,2,1,2), 2, 2, byrow=TRUE), c(2,2), respect=TRUE)
                            par(mar=c(4,4,4,0))
                            hdr.den(as.numeric(ALL.TFss.w2.arg.400x[m,3:402]), xlab=m, main="Histogram TF LOD m=ARG")
                            par(mar=c(4,2,4,0))
                            hdr.den(as.numeric(ALL.TFss.w2.imp.400x[m,3:402]), xlab=m, main="Histogram TF LOD m=IMP")
                            })
dev.off()


nf<-layout(matrix(c(1,2,1,2), 2, 2, byrow=TRUE), c(2,2), respect=TRUE)
layout.show(nf)
hdr.den(as.numeric(ALL.TFss.w2.arg.400x["D2Mit213",3:402]), xlab="LOD", main="HDR D2Mit213\n method=arg")
hdr.den(as.numeric(ALL.TFss.w2.imp.400x["D2Mit213",3:402]), xlab="LOD", main="HDR D2Mit213\n method=imp")

hdrs<- function(m){
      nf<-layout(matrix(c(1,2,1,2), 2, 2, byrow=TRUE), c(2,2), respect=TRUE)
      #layout.show(nf)   
        hdr.den(as.numeric(ALL.TFss.w2.arg.400x[m,3:402]), xlab="LOD", main=paste("HDR", m,"\n", "method=arg"))
        hdr.den(as.numeric(ALL.TFss.w2.imp.400x[m,3:402]), xlab="LOD", main=paste("HDR", m,"\n", "method=imp"))
                   }

pdf(file="../figures/Marker HDR density 400x.pdf")            
working.markers<- c(markers[1:151],markers[155:300],markers[301:450],markers[451:552], markers[554:557],markers[559:593])
sapply(working.markers, function(x) hdrs(x))
dev.off()

