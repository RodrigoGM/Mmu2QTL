###############################################################
# Replicated Composite Interval Mapping
#   Code in this file will be implemented in the composite
# interval mapping analysis. This analysis is described in
# further detail in Verdugo, RA (2007) THESIS.  Three step
# sizes c(1,.5, .25) and three window sizes c(1, .5, .25)
# will be used, in addition to the previous analysis with
# a step size = .1 and a window size of 2
#
#   Rodrigo Gularte-Mérida
#   9/24/2008 9:58:07 AM
###############################################################

setwd("../data")

require(qtl)|| stop('qtl package not available')

# read cross
D2sq <-  read.cross("csvs", genfile = "SBC_hg2d_geno.csv", phefile = "SBC_hg2d_pheno.csv" , na.strings = c(NA, ".", "-"), genotypes = c("B","H","C"), alleles = c("B","C"))

# new variable, lbwt: lean body weight
D2sq$pheno$lbwt <- D2sq$pheno$sac-D2sq$pheno$tf

# phenotype names
phenos<- c("sac", "gfp", "rfp", "mfp", "ffp", "tf", "lbwt")

# adjust fat traits for strain
# model: <phenotype> ~ strain 
phen.str <-  apply(D2sq$pheno[, phenos], 2,  function(p) {
  residuals(lm(p ~ D2sq$pheno$strain, na.action = na.exclude))
})
colnames(phen.str) <- paste(phenos, ".e", sep = "") # column names

# adjusting fat traits for strain, sac and sex to be used in CIM.
# model: <phenotype> ~ sex + sac
phen.ss <-  apply(D2sq$pheno[,phenos], 2,  function(p) {
  residuals(lm(p ~ D2sq$pheno$strain + D2sq$pheno$sex + D2sq$pheno$sac, na.action = na.exclude))
})
colnames(phen.ss) <- paste(phenos, ".ss", sep = "")  # column names

# adjusting fat traits for strain, sac, sex and sex*sac to be used in CIM.
# model: <phenotype> ~ sex + sac + sex:sac
phen.ss2 <-  apply(D2sq$pheno[,phenos], 2,  function(p) {
  residuals(lm(p ~ D2sq$pheno$strain + D2sq$pheno$sex + D2sq$pheno$sac + D2sq$pheno$sex:D2sq$pheno$sac,
               na.action = na.exclude))
})
colnames(phen.ss2) <- paste(phenos, ".ss2", sep = "")  # column names

# merging adjusted data to `<cross>$pheno`
D2sq$pheno <-  cbind(D2sq$pheno, phen.str, phen.ss, phen.ss2)

# estimating the recombination fraction of the genotype data
D2sq <-  est.rf(D2sq)
newmap <-  est.map(D2sq, error.prob = 0.01, verbose = TRUE)   # estimates a new gentic map

#estimates the genotype probabilities and a joint distribution
D2sq.4s <-  lapply(c(0, 1,.5,.25), function(step) {
                    calc.genoprob(D2sq, step = step, error.prob = 0.01)
                  })

# estimation of pseudo markers at different cM intervals
steps <- c(0, 1, 0.5, 0.25)

# split the simulation of pseudo markers due to memmory constraints
D2sq.4s. <-  lapply(1:2, function(cross.step) sim.geno(D2sq.4s[[cross.step]],
                    step = steps[cross.step], n.draws = 125, error.prob = 0.01))

D2sq.4s.. <-  lapply(3:4, function(cross.step) sim.geno(D2sq.4s[[cross.step]],
                    step = steps[cross.step], n.draws = 125, error.prob = 0.01))

( D2sq.s0 <-  D2sq.4s.[[1]] )
( D2sq.s1 <-  D2sq.4s.[[2]] )
( D2sq.s05 <-  D2sq.4s..[[1]] )
( D2sq.s025 <-  D2sq.4s..[[2]] )

rm(D2sq.4s., D2sq.4s..)

(crosses <- ls(pattern = "D2sq"))

# quality control
all.equal(D2sq.s0, D2sq.s1)
all.equal(D2sq.s0$pheno, D2sq.s1$pheno)
all.equal(D2sq.s0$geno, D2sq.s1$geno)
D2sq.s0$geno
D2sq.s1$geno
D2sq.s05$geno

D2sq.s025 <- calc.errorlod(D2sq.s025, error.prob = 0.01)
D2sq.s05 <- calc.errorlod(D2sq.s05, error.prob = 0.01)
D2sq.s1 <- calc.errorlod(D2sq.s1, error.prob = 0.01)
D2sq.s0 <- calc.errorlod(D2sq.s0, error.prob = 0.01)

# Saved each cross.  The s[0-9]+ indicate the Step in that cross
# an object named crosses was generated to hold the names of the crosses.
save(list = crosses, file = "D2sq_s0,s0.25,s0.5,s1.rda")  

