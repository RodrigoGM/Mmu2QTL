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
#   10/1/2008 4:47:04 PM
###############################################################

setwd("../analysis/")

require(qtl)|| stop('qtl package not available')

##DATA
# load indifidual cross files

TF.s0 <- matrix(NA, nrow=length(rownames(all.tf.s0.w1)), ncol=length(ls(pattern='.s0\\.')))
TF.s0 <- data.frame(TF.s0)
rownames(TF.s0) <- rownames(all.tf.s0.w1)
colnames(TF.s0) <- gsub('all\\.', '', ls(pattern='.s0\\.'))

TF.s1 <- matrix(NA, nrow=length(rownames(all.tf.s1.w1)), ncol=length(ls(pattern='.s1\\.')))
TF.s1 <- data.frame(TF.s1)
rownames(TF.s1) <- rownames(all.tf.s1.w1)
colnames(TF.s1) <- gsub('all\\.', '', ls(pattern='.s1\\.'))

TF.s05 <- matrix(NA, nrow=length(rownames(all.tf.s05.w1)), ncol=length(ls(pattern='tf.s05\\.')))
TF.s05 <- data.frame(TF.s05)
rownames(TF.s05) <- rownames(all.tf.s05.w1)
colnames(TF.s05) <- gsub("all\\.", "", ls(pattern='tf.s05\\.'))

sims <- colnames(all.tf.s0.w1)[grep('Sim', colnames(all.tf.s0.w1))]

## D2sq.s0
TF.s0[,'tf.s0.w1'] <- apply(all.tf.s0.w1[,sims], 1, median)
TF.s0[,'tf.s0.w05'] <- apply(all.tf.s0.w05[,sims], 1, median)
TF.s0[,'tf.s0.w025'] <- apply(all.tf.s0.w025[,sims], 1, median)

TF.s0 <- data.frame(chr=all.tf.s0.w1$chr, pos=all.tf.s0.w1$pos, TF.s0)

class(TF.s0) <- c('scanone', 'data.frame')
attr(TF.s0, "method") = 'hk'
#load("C:/.../D2sqs0")
attr(TF.s0, 'map') = pull.map(D2sq.s0)


## D2sq.s1
TF.s1[,'tf.s1.w1'] <- apply(all.tf.s1.w1[,sims], 1, median)
TF.s1[,'tf.s1.w05'] <- apply(all.tf.s1.w05[,sims], 1, median)
TF.s1[,'tf.s1.w025'] <- apply(all.tf.s1.w025[,sims], 1, median)

TF.s1 <- data.frame(chr=all.tf.s1.w1$chr, pos=all.tf.s1.w1$pos, TF.s1)

class(TF.s1) <- c('scanone', 'data.frame')
attr(TF.s1, "method") = 'hk'
#load("C:/.../D2sqs1")
attr(TF.s1, 'map') = pull.map(D2sq.s1)

## D2sq.s05
TF.s05[,'tf.s05.w1'] <- apply(all.tf.s05.w1[,sims], 1, median)
TF.s05[,'tf.s05.w05'] <- apply(all.tf.s05.w05[,sims], 1, median)
TF.s05[,'tf.s05.025'] <- apply(all.tf.s05.025[,sims], 1, median)

TF.s05 <- data.frame(chr=all.tf.s05.w1$chr, pos=all.tf.s05.w1$pos, TF.s05)

class(TF.s05) <- c('scanone', 'data.frame')
attr(TF.s05, "method") = 'hk'
#load("C:/.../D2sqs05")
attr(TF.s05, 'map') = pull.map(D2sq.s05)

legends <- names(TF.s0)[grep('s0', names(TF.s0))]
legends <- gsub('\\.', ' -', legends)

plot(TF.s05, lodcolum=c(1,2,3), xlim=c(145, 180), col=c(7,8,9))
plot(TF.s0, lodcolum=c(1,2,3), xlim=c(145, 180), col=c(1,2,3), add=TRUE)
plot(TF.s1, lodcolum=c(1,2,3), xlim=c(145, 180), col=c(4,5,6), add=TRUE)

plot(TF.s05, lodcolum=c(1,2,3), xlim=c(172, 176), col=c("purple"), lwd=2, lty=c(1,2,3), ylim=c(-1, 3.5), incl.markers=TRUE, show.marker.names=TRUE)
plot(TF.s0, lodcolum=c(1,2,3), xlim=c(173, 175), col=c("brown"), lwd=2, lty=c(1,2,3), add=TRUE)
plot(TF.s1, lodcolum=c(1,2,3), xlim=c(173, 175), col=c("darkgreen"),lwd=2, lty=c(1,2,3), add=TRUE)


