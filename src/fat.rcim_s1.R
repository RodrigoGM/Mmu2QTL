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
#   8.Oct.2008
###############################################################
# libraries
require(qtl)
source("../src/qtl.R")
source("../src/qtl2.R")

# working directory
setwd('../analysis/')

# data
load("../data/D2sq_s0,s0.25,s0.5,s1.rda")

# code profiling
Rprof()

# set up window sizes
winsize <- c(1, .5, .25)

# run rCIM for tota fat mass
rcim.tf.ss.s1 <- lapply(winsize, function(win.size) {
	r.cim(D2sq.s1, pheno.col = c('tf.ss'), times = 400, n.marcovar = 3, windo = win.size, method = 'hk', save = TRUE)

	if(win.size == 1) {all.tf.s1.w1 <<- rcim.results}
	if(win.size == .5) {all.tf.s1.w05 <<- rcim.results}
	if(win.size == .25) {all.tf.s1.w025 <<- rcim.results}

	})

#save and discard the results after the loop
save(list = c(ls(pattern = 'all.tf'),'rcim.tf.ss.s1'), file = "rcim/tf_s1_w1_w0.5_w0.25.rda")
rm(list = c(ls(pattern = 'all.tf'),'rcim.tf.ss.s1'))

# run rCIM for independent fat pad tissues
# save and discard results at the end of each loop
rcim.gfp.ss.s1 <- lapply(winsize, function(win.size) {
	r.cim(D2sq.s1, pheno.col = c('gfp.ss'), times = 400, n.marcovar = 3, window = win.size, method = 'hk', save = TRUE)

	if(win.size == 1) {all.gfp.s1.w1 <<- rcim.results}
	if(win.size == .5) {all.gfp.s1.w05 <<- rcim.results}
	if(win.size == .25) {all.gfp.s1.w025 <<- rcim.results}

	})

save(list = c(ls(pattern = 'all.gfp'),'rcim.gfp.ss.s1'), file = "rcim/gfp_s1_w1_w0.5_w0.25.rda")
rm(list = c(ls(pattern = 'all.gfp'),'rcim.gfp.ss.s1'))

rcim.rfp.ss.s1 <- lapply(winsize, function(win.size) {
	r.cim(D2sq.s1, pheno.col = c('rfp.ss'), times = 400, n.marcovar = 3, window = win.size, method = 'hk', save = TRUE)

	if(win.size == 1) {all.rfp.s1.w1 <<- rcim.results}
	if(win.size == .5) {all.rfp.s1.w05 <<- rcim.results}
	if(win.size == .25) {all.rfp.s1.w025 <<- rcim.results}

	})

save(list = c(ls(pattern = 'all.rfp'),'rcim.rfp.ss.s1'), file = "rcim/rfp_s1_w1_w0.5_w0.25.rda")
rm(list = c(ls(pattern = 'all.rfp'),'rcim.rfp.ss.s1'))


rcim.ffp.ss.s1 <- lapply(winsize, function(win.size) {
	r.cim(D2sq.s1, pheno.col = ('ffp.ss'), times = 400, n.marcovar = 3, window = win.size, method = 'hk', save = TRUE)

	if(win.size == 1) {all.ffp.s1.w1 <<- rcim.results}
	if(win.size == .5) {all.ffp.s1.w05 <<- rcim.results}
	if(win.size == .25) {all.ffp.s1.w025 <<- rcim.results}

	})

save(list = c(ls(pattern = 'all.ffp'),'rcim.ffp.ss.s1'), file = "rcim/ffp_s1_w1_w0.5_w0.25.rda")
rm(list = c(ls(pattern = 'all.ffp'),'rcim.ffp.ss.s1'))


rcim.mfp.ss.s1 <- lapply(winsize, function(win.size) {
	r.cim(D2sq.s1, pheno.col = ('mfp.ss'), times = 400, n.marcovar = 3, window = win.size, method = 'hk', save = TRUE)

	if(win.size == 1) {all.mfp.s1.w1 <<- rcim.results}
	if(win.size == .5) {all.mfp.s1.w05 <<- rcim.results}
	if(win.size == .25) {all.mfp.s1.w025 <<- rcim.results}

	})

save(list = c(ls(pattern = 'all.mfp'),'rcim.mfp.ss.s1'), file = "rcim/mfp_s1_w1_w0.5_w0.25.rda")
rm(list = c(ls(pattern = 'all.mfp'),'rcim.mfp.ss.s1'))

