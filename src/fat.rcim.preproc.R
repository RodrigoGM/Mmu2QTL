###############################################################
## Replicated Composite Interval Mapping                     ##
##   Code in this file will be implemented in the composite  ##
## interval mapping analysis. This analysis is described in  ##
## further detail in Verdugo, RA (2007) THESIS.  Three step  ##
## sizes c(1,.5, .25) and three window sizes c(1, .5, .25)   ##
## will be used, in addition to the previous analysis with   ##
## a step size = .1 and a window size of 2  -R Gularte       ##
###############################################################

setwd("../../data")

################################
##    9/24/2008 9:58:07 AM    ##
################################

require(qtl)|| stop('qtl package not available')

# READ CROSS
D2sq<- read.cross("csv", dir="C:\\Documents and Settings\\RJGularte\\My Documents\\LAB FILES\\2D Congenics\\2D Sub-congenics\\QTL\\CIM",
        file=list.files(pattern="FINAL SET\\.csv"),                    #  Base File
        na.strings=c(".","NA"), genotypes=c("B","H","C"), alleles=c("B","C"))

D2sq$pheno<- cbind(D2sq$pheno, lbwt=D2sq$pheno$SAC-D2sq$pheno$TF)                                 #new variable, lbwt: lean body weight

# NAMES OF THE PHENOTYPES
phenotypes<-c(names(D2sq$pheno[find.pheno(D2sq, "WK2"):find.pheno(D2sq, "Muscle")]),'lbwt')

# CORRECT DATA BY STRAIN ONLY
phen.str<- data.frame(apply(D2sq$pheno[, phenotypes], 2,  function(p) {
                  residuals(lm(p ~ D2sq$pheno$Line,                                                                    #model: <phenotype> ~ strain
                  na.action=na.exclude))}))                                                                            #data.frame of residuals only by strain
colnames(phen.str)<-paste(phenotypes, ".e", sep="")                                                                    #set column names
#D2sq$pheno<- cbind(D2sq$pheno, phen.res)                                                                              # New 'pheno' data.frame with corrected data

## Adjusting data for Strain, sac and sex to be used in CIM.
phen.ss<- data.frame(apply(D2sq$pheno[,phenotypes], 2,  function(p) {
                residuals(lm(p ~ D2sq$pheno$Line + D2sq$pheno$SEX + D2sq$pheno$SAC,      #model: <phenotype> ~ sex + sac
                na.action=na.exclude))}))
colnames(phen.ss)<-paste(phenotypes, ".ss", sep="")  #colum names
#D2sq$pheno<- cbind(D2sq$pheno, phen.ss)

## Adjusting data for Strain, sac, sex and sex*sacto be used in CIM.
phen.ss2<- data.frame(apply(D2sq$pheno[,phenotypes], 2,  function(p) {
                  residuals(lm(p ~ D2sq$pheno$Line + D2sq$pheno$SEX + D2sq$pheno$SAC + D2sq$pheno$SEX:D2sq$pheno$SAC,
                  na.action=na.exclude))}))
colnames(phen.ss2)<-paste(phenotypes, ".ss2", sep="")  #colum names
#phen.ss<-data.frame(phen.res) # corrected data
#D2sq$pheno<- cbind(D2sq$pheno, phen.ss2)

## Adjusting WK6, WK9, ECW and lbwt for Strain, wk2 and sex to be used in CIM
phen.sw<- data.frame(apply(D2sq$pheno[,c('WK6', 'WK9', 'ECW', 'lbwt')], 2,  function(p) {
                  residuals(lm(p ~ D2sq$pheno$Line + D2sq$pheno$SEX + D2sq$pheno$WK2,
                  na.action=na.exclude))}))
colnames(phen.sw)<-c("WK6.sw", "WK9.sw", "ECW.sw", "lbwt.sw")           #colum names
#D2sq$pheno<- cbind(D2sq$pheno, phen.sw)

##  Adjusting WK6, WK9, ECW and lbwt for Strain, wk2, sex and age
phen.a<- data.frame(apply(D2sq$pheno[,c('WK6', 'WK9', 'ECW', 'lbwt')], 2,  function(p) {
                  residuals(lm(p ~ D2sq$pheno$Line + D2sq$pheno$SEX + D2sq$pheno$WK2 + D2sq$pheno$Age,
                  na.action=na.exclude))}))
colnames(phen.a)<-c("WK6.a", "WK9.a", "ECW.a", "lbwt.a")           #colum names

D2sq$pheno<- cbind(D2sq$pheno, phen.str, phen.ss, phen.ss2, phen.sw, phen.a)

# ESTIMATING THE RECOMBINATION FRACTION
D2sq<- est.rf(D2sq)         #estimating the recombination fraction of the data
newmap<- est.map(D2sq, error.prob=0.01, verbose=TRUE)   # estimates a new gentic map

save.image()
# ESTIMATE GENOTYPE PROBABILITIES AND SIMULATE PSEUDOMARKERS
D2sq.4s<- lapply(c(0, 1,.5,.25), function(step) {
                    calc.genoprob(D2sq, step=step, error.prob=0.01)
                 }) #estimates the genotype probabilities and a joint distribution
save.image()

steps<-c(0, 1,.5,.25)
D2sq.4s.<- lapply(1:2, function(cross.step) sim.geno(D2sq.4s[[cross.step]],
                    step=steps[cross.step], n.draws=125, error.prob=0.01))

D2sq.4s..<- lapply(3:4, function(cross.step) sim.geno(D2sq.4s[[cross.step]],
                    step=steps[cross.step], n.draws=125, error.prob=0.01))

D2sq.s0<- D2sq.4s.[[1]]
D2sq.s1<- D2sq.4s.[[2]]
D2sq.s05<- D2sq.4s..[[1]]
D2sq.s025<- D2sq.4s..[[2]]
rm(D2sq.4s., D2sq.4s..)
ls()
crosses<-ls()
crosses
all.equal(D2sq.s0, D2sq.s1)
all.equal(D2sq.s0$pheno, D2sq.s1$pheno)
all.equal(D2sq.s0$geno, D2sq.s1$geno)
D2sq.s0$geno
D2sq.s1$geno
D2sq.s05$geno
ls()


# Saved each cross as an object.  The s[0-9]+ indicate the Step in that cross
# an object named crosses was generated to hold the names of the crosses.
### save(list=ls(), file="08-09-25 D2sq steps 0 1 05 025, 4 obj.RData")  


lapply(c(D2sq.s0, D2sq.s1, D2sq.s05, D2sq.s025), function(crossi) crossi<<-calc.errorlod(crossi, error.prob=0.01))   # Calculating errors in genotyping

D2sq.s025<-calc.errorlod(D2sq.s025, error.prob=0.01)
D2sq.s05<-calc.errorlod(D2sq.s05, error.prob=0.01)
D2sq.s1<-calc.errorlod(D2sq.s1, error.prob=0.01)
D2sq.s0<-calc.errorlod(D2sq.s0, error.prob=0.01)

winsize<-c(1, .5, .25) 

