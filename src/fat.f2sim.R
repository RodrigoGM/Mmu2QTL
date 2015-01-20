#################################################################
# CROSS SIMULATION USING THE ESTIMATED POSITIONS AND USING 
#  ADDITIVE AND DOMINANCE EFFECTS OF MARKERS RS... AND D2MIT213 
#  FROM ALL THE SUBCONGENIC STRAINS
#  SASglm: TF = A1 D1 SEX SAC LINE; TF =  A2 D2 SEX SAC LINE;
#  SIMULATION OF CROSS END.
# DETAILS
#     The sim.1 is similar to the genotype information provided by a microsatellite panel
#     The sim.2 is similar to the genptype information provided by a combination of SNP and Microsatellite panels
#     Each sim.# has the same phenotypic information, and a control; values obtained from a normal distribution 
#       with the same mean, variance, additive and dominance values of the Fatq2 QTL, but with no specific location 
#     The effects of both QTL were given by the a and d of Total Fat from the HG2D- subcongenics panel where Q1 is the most distal QTL
#       and Q2 is the minor QTL located around Agouti.
#################################################################
## WORKING DIRECTORY
setwd("../analysis/")
## LIBRARIES
library(qtl)
library(coda)
library(hdrcde)
library(boa)
source('../src/qtl.R')
source('../src/qtl2.R')


## DATA AND OUTPUTS
load("../data/SBC_hg2d.rda")
load("scanone_output.rda")
#load("rcim simulation results.rda")

#load(file="08-04-29 sim.1.rda")      #  Simulation 1: 11 markers not equally spaced e.g. microsattelites
#load(file="08-04-29 sim.2.rda")      #  Simulation 2: 53 markers equally spaced e.g. SNP
#load("08-04-30 CIM of control phenotype.rda")  # Simulations 1 and 2
#load("08-04-30 rCIM of Alternative phenotype.rda")


# seed for random number generator
seed <- ceiling(runif(1, 0, 10^8))
seed
  #[1] 4914561
set.seed(seed)

# simulating maps
map <- sim.map(len = 50, n.mar = 11, include.x = FALSE, eq.spacing = TRUE)
map2 <- sim.map(len = 50, n.mar = (11+42), include.x = FALSE, eq.spacing = TRUE)
par(mfrow=c(1,2))
plot(map); plot(map2)

# simulate cross
sim.1 <- sim.cross(map, type = "f2", n.ind = 1920, model = rbind(c(1,42.1,-.8,-.004),c(1,20,-0.05,-0.007)), error.prob = .001, missing.prob = 0.009)
sim.2 <- sim.cross(map2, type = "f2", n.ind = 1920, model = rbind(c(1,42.1,-.8,-.004),c(1,20,-0.05,-0.007)), error.prob =  .001, missing.prob = 0.009)

all.equal(sim.1, sim.2)
cbind(sim.1$pheno, sim.2$pheno)

# Adapted from Siegmund and Yakir 2007; page 41
n <- 1920; p <-0.05; x <- rbinom(n, 3, p)
mu1 <- 1.097865; a2 <- -0.08; d2 <- -0.004; sig1 <- .192
mu2 <- 1.097865; a1 <- -0.05; d1 <- -0.007; sig2 <- .194

y <- mu1 + a1*x + d1*(x==1) + rnorm(n, sd=sig1)
y2 <- mu2 + a2*x + d2*(x==1) + rnorm(n, sd=sig2)
Ctrl <- (y+y2)/2

sim.1$pheno <- cbind(sim.1$pheno, Ctrl, pheno.a =  sim.2$pheno$phenotype)
sim.2$pheno <- cbind(sim.2$pheno, Ctrl, pheno.a =  sim.1$pheno$phenotype)

# Estimation of the recombination frequency and genetic map
sim.1 <- est.rf(sim.1)
sim.2 <- est.rf(sim.2)
newmap.1 <- est.map(sim.1, error.prob = 0.01, verbose = TRUE)
newmap.2 <- est.map(sim.2, error.prob = 0.01, verbose = TRUE)

# calculation of genotype probabilities, genotyping errors, and pseudomarkers
# sim.1 <- calc.genoprob(sim.1, error.prob = 0.0001)
# sim.2 <- calc.genoprob(sim.2, error.prob = 0.0001)
sim.1 <- calc.genoprob(sim.1, step = .1, error.prob = 0.01) 
sim.2 <- calc.genoprob(sim.2, step = .1, error.prob = 0.01) 
sim.1 <- sim.geno(sim.1, n.draws = 100, step = 1, error.prob = 0.0001)
sim.2 <- sim.geno(sim.2, n.draws = 100, step = .5, error.prob = 0.0001)
sim.1 <- argmax.geno(sim.1, error.prob = 0.01)
sim.2 <- argmax.geno(sim.2, error.prob = 0.01)
sim.1 <- calc.errorlod(sim.1, error.prob = 0.01)
sim.2 <- calc.errorlod(sim.2, error.prob = 0.01)


save(sim.1, file = "sim.1.rda")      # Saving Simulation 1: 11 markers not equally spaced
save(sim.2, file = "sim.2.rda")      # Saving Simulation 2: 53 markers equally spaced

#phenotypes <- names(sim.1$pheno)
out.scan.sim1 <- scanone(sim.1, pheno.col = c("phenotype", "Ctrl", "pheno.a"), method = "hk")
out.scan.sim2 <- scanone(sim.2, pheno.col = c("phenotype", "Ctrl", "pheno.a"), method = "hk")


## COMPOSITE INTERVAL MAPPING
## real phenotype (rp)
#set.seed(seed)
# sim1
sim1.w2.rp <- r.cim(sim.1, pheno.col = "phenotype", times = 300, method = "hk", save = TRUE)
all.sim1.w2.rp <-rcim.results
sim1.w1.rp <- r.cim(sim.1, pheno.col = "phenotype", times = 300, method = "hk", save = TRUE, window = 1)
all.sim1.w1.rp <- rcim.results
# sim2
sim2.w2.rp <- r.cim(sim.2, pheno.col = "phenotype", times = 300, method = "hk", save = TRUE)
all.sim2.w2.rp <-rcim.results
sim2.w1.rp <- r.cim(sim.2, pheno.col = "phenotype", times = 300, method = "hk", save = TRUE, window = 1)
all.sim2.w1.rp <- rcim.results


## control phenotype
#set.seed(seed)
sim1.w2.ct <- r.cim(sim.1, pheno.col = "Ctrl", times = 300, method = "hk", save = TRUE)
all.sim1.w2.ct <- rcim.results
sim1.w1.ct <- r.cim(sim.1, pheno.col = "Ctrl", times = 300, method = "hk", save = TRUE, window = 1)
all.sim1.w1.ct <- rcim.results
sim2.w2.ct <- r.cim(sim.2, pheno.col = "Ctrl", times = 300, method = "hk", save = TRUE)
all.sim2.w2.ct <- rcim.results
sim2.w1.ct <- r.cim(sim.2, pheno.col = "Ctrl", times = 300, method = "hk", save = TRUE, window = 1)
all.sim2.w1.ct <- rcim.results

##save(sim1.w2, sim1.w1, all.sim1.w2, all.sim1.w1, sim2.w2, all.sim2.w2, sim2.w1, all.sim2.w1, file = "simulated_crosses_1+2_s_01_m_3.rda")
save(sim1.w1.ct, all.sim1.w1.ct, sim2.w1.ct, all.sim2.w1.ct, file = "control_phenotype.rda")


## analysis of alternate phenotype
#set.seed(seed)
sim1.w1.A <- r.cim(sim.1, pheno.col = "pheno.a", times = 300, method = "hk", save = TRUE, window = 1)
all.sim1.w1.A <- rcim.results
sim2.w1.A <- r.cim(sim.2, pheno.col = "pheno.a", times = 300, method = "hk", save = TRUE, window = 1)
all.sim2.w1.A <- rcim.results
sim1.w2.A <- r.cim(sim.1, pheno.col = "pheno.a", times = 300, method = "hk", save = TRUE)
all.sim1.w2.A <- rcim.results
sim2.w2.A <- r.cim(sim.2, pheno.col = "pheno.a", times = 300, method = "hk", save = TRUE)
all.sim2.w2.A <- rcim.results

sim2.w2.10m.rt <- r.cim(sim.2, pheno.col = "phenotype", n.marcovar = 10, times = 300, method = "hk", save = TRUE)
all.sim2.w2.10m.rt <- rcim.results

save(sim1.w1.A, all.sim1.w1.A, sim2.w1.A, all.sim2.w1.A, file = "alternate_phenotype.rda")

#save(sim1.w2, sim1.w1, all.sim1.w2, all.sim1.w1, sim2.w2, all.sim2.w2, sim2.w1, all.sim2.w1, file = "08-04-30 Simulated crosses W_1&2, S_01, M_3.rda")
#save(list = ls(), file = "08-04-30 Simulations w_ 11 and 53 markers, s_0.1 w_1 w_2.rda")


pdf("SimulatedLinkedQTL_rCIM.pdf", width = 297/25.4, height = 210/25.4)
par(mfrow = c(1,2))
plot(out.scan.sim1, lodcolumn = 1:3, col = c("red","orange","purple"), 
          main = "Simulated cross 1  w/QTLs at 20.1 and 42.1 cM\n 11 markers, equally spaced, step = .1, 300x CIM")
plot(sim1.w1.rp, sim1.w1.ct, sim1.w1.A, col = c("black","blue","brown"), add = TRUE)
plot(sim1.w2.rp, sim1.w2.ct, sim1.w2.A, col = c("pink","gray","magenta"), add = TRUE)
legend("topleft", c("HK-Scan: Phenotype","HK-Scan CTRL","HK-Scan Alternate",
      "rCIM: W = 1 Phenotype ", "rCIM: W = 1 CTRL ", "rCIM: W = 1 Alternate",
      "rCIM: W = 2 Phenotype ", "rCIM: W = 2 CTRL ", "rCIM: W = 2 Alternate"),
       col = c("red","orange", "purple", "black","blue","brown", "pink", "gray","magenta"), lty = 1, lwd = 2, bty = "n")

plot(out.scan.sim2, lodcolumn = 1:3, col = c("red","orange","purple"), 
          main = "Simulated cross 2 w/QTLs at 20.1 and 42.1 cM\n 53 markers, equally spaced, step = .1, 300x CIM")
plot(sim2.w1.rp, sim2.w1.ct, sim2.w1.A, col = c("blue", "black","brown"),  add = TRUE)
plot(sim2.w2.rp, sim2.w2.ct, sim2.w2.A, col = c("pink", "gray","magenta"), add = TRUE)
legend("topleft", c("HK-Scan Phenotype","HK-Scan CTRL","HK-Scan Alternate",
      "rCIM: W = 1 Phenotype", "rCIM: W = 1 CTRL", "rCIM W = 1 Alternate",
      "rCIM: W = 2 Phenotype", "rCIM: W = 2 CTRL", "rCIM W = 2 Alternative"),
       col = c("red","orange", "purple", "black","blue","brown", "pink", "gray","magenta"), lty = 1, lwd = 2, bty = "n")
dev.off()

plot(Sim.2 ~ pos, data = all.sim2.w2.rp, type = "l", ylim = c(0, 100))
invisible(
sapply(4:302, function(i) lines(x = all.sim2.w2.rp$pos, y = all.sim2.w2.rp[, i], type = "l"))
)

plot(Sim.2 ~ pos data = all.sim2.w2.10m.rt, type = 'l')
invisible(
sapply(4:302, function(i) lines(x = all.sim2.w2.10m.rt$pos, y = all.sim2.w2.10m.rt[, i], type = "l"))
)


cim.sim2.w2.10m.rt <- cim(sim.2, pheno.col = 1, n.marcovar = 10, method = "hk", window = 2)
cim.sim2.w10.10m.rt <- cim(sim.2, pheno.col = 1, n.marcovar = 10, method = "hk", window = 10)
cim.sim2.w2.3m.rt <- cim(sim.2, pheno.col = 1, n.marcovar = 3, method = "hk", window = 2)
cim.sim2.w10.3m.rt <- cim(sim.2, pheno.col = 1, n.marcovar = 3, method = "hk", window = 10)

plot(cim.sim2.w2.10m.rt, cim.sim2.w10.10m.rt, cim.sim2.w2.3m.rt, ylab = "LOD", xlab = "Map Position (cM)",
   main = "Simulated Cross 1 w/QTLs at 20.1 and 42.1 cM\n 53 markers, equally spaced, step = .1, 1x CIM")
plot(cim.sim2.w10.3m.rt,col = "green3",  add = TRUE)
legend("topleft", legend = c("CIM : W = 2, M = 10", "CIM : W = 10, M = 10", "CIM : W = 2, M = 3", "CIM : W = 10, M = 3"),  
     col = c("black", "blue", "red", "green3"), lty = 1)

