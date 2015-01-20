###############################################################
## Replicated Composite Interval Mapping                     ##
##   Code in this file will be implemented in the composite  ##
## interval mapping analysis. This analysis is described in  ##
## further detail in Verdugo, RA (2007) THESIS.  Three step  ##
## sizes c(1,.5, .25) and three window sizes c(1, .5, .25)   ##
## will be used, in addition to the previous analysis with   ##
## a step size = .1 and a window size of 2  -R Gularte       ##
###############################################################

setwd('/Server/Meyer_Transfer_Folder/Rodrigo/D2sqs025')

Rprof()

winsize<-c(1, .5, .25)

rcim.tf.ss<-lapply(winsize, function(win.size) {
	r.cim(D2sq.s025, pheno.col=c('TF.ss'), times=400, n.marcovar=3, window=win.size, method='hk', save=TRUE)

	if(win.size==1) {all.tf.s025.w1<<-rcim.results}
	if(win.size==.5) {all.tf.s025.w05<<- rcim.results}
	if(win.size==.25) {all.tf.s025.w025<<-rcim.results}

	})

save(list=c(ls(pattern='all.tf'),'rcim.tf.ss'), file="08-10-08 TF -s025 -w1,05,025.RData")
rm(list=c(ls(pattern='all.tf'),'rcim.tf.ss'))

rcim.gfp.ss<-lapply(winsize, function(win.size) {
	r.cim(D2sq.s025, pheno.col=c('GFP.ss'), times=400, n.marcovar=3, window=win.size, method='hk', save=TRUE)

	if(win.size==1) {all.gfp.s025.w1<<-rcim.results}
	if(win.size==.5) {all.gfp.s025.w05<<- rcim.results}
	if(win.size==.25) {all.gfp.s025.w025<<-rcim.results}

	})

save(list=c(ls(pattern='all.gfp'),'rcim.gfp.ss'), file="08-10-08 GFP -s025 -w1,05,025.RData")
rm(list=c(ls(pattern='all.gfp'),'rcim.gfp.ss'))

rcim.RFP.ss.sw<-lapply(winsize, function(win.size) {
	r.cim(D2sq.s025, pheno.col=c('RFP.ss'), times=400, n.marcovar=3, window=win.size, method='hk', save=TRUE)

	if(win.size==1) {all.rfp.s025.w1<<-rcim.results}
	if(win.size==.5) {all.rfp.s025.w05<<- rcim.results}
	if(win.size==.25) {all.rfp.s025.w025<<-rcim.results}

	})

save(list=c(ls(pattern='all.rfp'),'rcim.rfp.sw'), file="08-10-08 RFP -s025 -w1,05,025.RData")
rm(list=c(ls(pattern='all.rfp'),'rcim.rfp.sw'))


rcim.ffp.ss<-lapply(winsize, function(win.size) {
	r.cim(D2sq.s025, pheno.col=c('FFP.ss'), times=400, n.marcovar=3, window=win.size, method='hk', save=TRUE)

	if(win.size==1) {all.ffp.s025.w1<<-rcim.results}
	if(win.size==.5) {all.ffp.s025.w05<<- rcim.results}
	if(win.size==.25) {all.ffp.s025.w025<<-rcim.results}

	})

save(list=c(ls(pattern='all.ffp'),'rcim.ffp.ss'), file="08-10-08 FFP -s025 -w1,05,025.RData")
rm(list=c(ls(pattern='all.ffp'),'rcim.ffp.ss'))


rcim.mfp.ss<-lapply(winsize, function(win.size) {
	r.cim(D2sq.s025, pheno.col=c('MFP.ss'), times=400, n.marcovar=3, window=win.size, method='hk', save=TRUE)

	if(win.size==1) {all.mfp.s025.w1<<-rcim.results}
	if(win.size==.5) {all.mfp.s025.w05<<- rcim.results}
	if(win.size==.25) {all.mfp.s025.w025<<-rcim.results}

	})

save(list=c(ls(pattern='all.ffp'),'rcim.ffp.ss'), file="08-10-08 MFP -s02 -w1,05,025.RData")
rm(list=c(ls(pattern='all.ffp'),'rcim.ffp.ss'))


save.image()
