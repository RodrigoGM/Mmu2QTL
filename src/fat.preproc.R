#################################################################
# The code written in this file is developed for the fine
#  mapping analysis of Fatq2 using SNP and microsattelite
#  as markers
#  The SBC_hg2d.rda data set contains the phenotype and genotype
#  from all five subcongenics and the data from the HG2D cross.
#  All genotypes and phenotypes were manually cuarated.  Three
#  SNP were removed because they lacked quality peaks (calls).
# To cite this data please use doi:10.5281/zenodo.12793, and the
#  publication.
#################################################################

# set data as working directory
setwd("../data/")
# require qtl
require(qtl) || stop("R/qtl is not available")


# reading data
D2sq<- read.cross("csvs", genfile = "SBC_hg2d_geno.csv", phefile = "SBC_hg2d_pheno.csv" , na.strings=c(NA, ".", "-"), genotypes=c("B","H","C"), alleles=c("B","C"))

# phenotypes to analyze
phenos<- c("sac", "gfp", "rfp", "mfp", "ffp", "tf")

# Adjusting phenotypes for strain effects, and binding the data to phenotype table in the cross object
phen.res<- data.frame(apply(D2sq$pheno[ ,phenos], 2,  function(p) {residuals(lm(p ~ D2sq$pheno$strain, na.action=na.exclude))}))
colnames(phen.res)<-paste(phenos, ".e", sep="")
D2sq$pheno<- cbind(D2sq$pheno, phen.res)

# Adjusting data for strain, sacrifice weight and sex
phen.ss<- data.frame(apply(D2sq$pheno[, phenos], 2,  function(p) {residuals(lm(p ~ D2sq$pheno$strain + D2sq$pheno$sex + D2sq$pheno$sac, na.action=na.exclude))}))
colnames(phen.ss)<-paste(phenos, ".ss", sep="")
D2sq$pheno<- cbind(D2sq$pheno, phen.ss)

# Adjusting data for strain, sacrifice weight, sex and sex*sac
phen.ss2<- data.frame(apply(D2sq$pheno[, phenos], 2,  function(p) {residuals(lm(p ~ D2sq$pheno$strain + D2sq$pheno$sex + D2sq$pheno$sac + D2sq$pheno$sex:D2sq$pheno$sac, na.action=na.exclude))}))
colnames(phen.ss2)<-paste(phenos, ".ss2", sep="")
D2sq$pheno<- cbind(D2sq$pheno, phen.ss2)


# Estimating the recombination fraction of the genotypes
D2sq<- est.rf(D2sq)

# Estimating genotype probabilities and simulating pseudomarkers
D2sq<- calc.genoprob(D2sq, step=.2, error.prob=0.01)
D2sq<- calc.errorlod(D2sq, error.prob=0.01)
D2sq<- sim.geno(D2sq, step=0, n.draws=125, error.prob=0.01)
D2sq<- argmax.geno(D2sq, error.prob=0.01)

# saves D2sq and phenotypes for later use.  Files will not be read, rather the rda will be loaded
save(D2sq, phenos, file="SBC_hg2d.rda")
