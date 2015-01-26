#################################################################
# Figure 2.pdf, Figure 3.pdf, Figure 4.pdf
# 
# Figure 2. Plots the fat QTL, identified in this manuscript, and
#  the locations of those previously identified by QTL from
#  Farber 2007b that are relevant to this work. We modified the
#  pdf generated with this script using Inkscape, modifications
#  can be seen in the corresponding SVG for Figure 2.
#################################################################

setwd("../figures")
library(qtl)
library(RColorBrewer)
source("../src/qtl.R")
source("../src/qtl2.R")

# load outputs and data
load("../analysis/scanone_output.rda")
load("../analysis/cim_output.rda")
load("../analysis/TF_rcimSummary.rda")
load("../data/SBC_hg2d.rda")

# additional data
(abv <- read.delim("../data/abv.dat", row.names = 1))
(markers <- names(pull.map(D2sq)[[1]]))
(qtls<- read.csv("../data/Chr2QTL_CRF.csv", row.names = 1))
(congenics <- read.delim("../data/congenic_strains.txt"))

# Figure 2
# 9/26/2008 11:53:18 AM

pdf("../figures/Figure 2Test.pdf", height = 210/25.4, width = 297/25.4)

nf <- layout(mat=matrix(c(1,2,1,3,4,3), nrow=3, ncol=2, byrow=TRUE),
             widths=c(3,1), heights=c(1,1,.7), respect=TRUE)

cols = brewer.pal(5, "Set1")[c(1,2,3,5,4)]

par(mar=c(0, 5, 3, 3)+.1, bty='l', xaxt='s', yaxt='s', cex.lab=2, cex.axis=1.8, lwd = 2)

plot.mscanone(out.FATtrte.ss, xlim=c(75,180), lwd=4, ylab="LOD", xlab="Map Position (Mbp)",
              col=cols, mtitle=FALSE,
              title="", incl.legend = FALSE, incl.markers = TRUE)
abline(h=3.5, col='gray', lwd=2, lty = 2)
text(x = 78, y = 3.5, labels = expression(alpha~"="~0.01), pos = 3, cex = 1.7, font = 2)
abline(h=2.5, col="gray1", lwd=2, lty = 2)
text(x = 78, y = 2.5, labels = expression(alpha~"="~0.05), pos = 3, cex = 1.7, font = 2)
mtext("A", side = 3, line = 1, outer = FALSE, at = 70, font = 2, cex = 1.8)

text(x = 168, y = 13, labels="Fatq2b", cex=1.9, pos = 2, font = 3)
text(x = 150, y = 6, labels = "Fatq2a", cex = 1.9, pos = 2, font = 3)
lines(x = c(152.5, 159), y = c(6.0, 6.0)-.3, col = cols[5], lwd = 3)
lines(x = c(151.9, 160), y = c(6.0, 6.0)+0, col = cols[1], lwd = 3)
lines(x = c(150, 161), y = c(6.0, 6.0)+.3, col = cols[2], lwd = 3)
lines(x = c(150, 161), y = c(6.0, 6.0)+.6, col = cols[3], lwd = 3)

legend('topleft',  c('FFP','TF', 'GFP', 'RFP', 'MFP'),
       col=cols[c(5,1,2,3,4)],
       lty=c(1,1,1,1), lwd=4.25, bty='n', cex=1.7)

# effect plots
par(mar=c(1, 5.9, 3, 0)+.1, xaxt='s', yaxt='s', bty='l', cex.lab=1.7, cex.axis=1.7, lwd = 2)
rs27669516<- find.marker(D2sq, chr=2, pos=174)
effectplot(D2sq, pheno.col=find.pheno(D2sq, "tf.ss"), mname1=rs27669516, main="", ylab="Total Fat (g)\n")
text(x=1, y=-.08, labels="174.0 Mb", cex=1.9, adj=c(.1,0), font = 2)
mtext("B", side = 3, line = 1, outer = FALSE, at = 0, font = 2, cex = 1.8)

par(mar=c(13.5, 5.9, 5.5, 0)+.1, xaxt='s', yaxt='s', bty='l', cex.lab=1.6, cex.axis=1.7, lwd = 2)
rs27307998<- find.marker(D2sq, chr=2, pos=156.97617)
effectplot(D2sq, pheno.col=find.pheno(D2sq, "tf.ss"), mname1=rs27307998, main="", ylab="Total Fat (g)\n", xlab="Genotype")
text(x=1, y=-.062, labels="156.9 Mb", cex=1.9, adj=c(.1,0), font = 2)
mtext("C", side = 3, line = 1, outer = FALSE, at = 0, font = 2, cex = 1.8)


# plot congenic strains on panel 5

par(mar=c(5, 5, 0, 3)+.1, yaxt='n', xaxt='s', cex.lab=1.7, cex.axis=1.6)
plot(out.TF.nc$lod ~ out.TF.nc$pos, type='n', ylim=c(-.13,1.2), xlim=c(75,180),
     xlab="Map Position (Mbp)", ylab="", bty='n')
mtext("D", side = 2, line = 0, cex = 1.8, font = 2, las = 2, padj = 1)

# coordinates to plot
congenics<- cbind(congenics, y=c(1/3,2/3,2/3,1/3,2/3, NA))

lapply(1:5, function(i) {
  lines(x=c(congenics[i,'start'],congenics[i,'end']), y=c(congenics[i,'y'],congenics[i,'y']),  # coordinates
        type='l', pch="", lwd=4, cex=1, col='black')            # line parameters
  text(x=congenics[i,'start'], y=congenics[i,'y']-.01,        #coordinates
       labels=congenics[i,'congenic'], pos = 1, cex=1.6) #, adj=c(0,1.1))
})

# Fatq1
lines(x=c(qtls['Fatq1', 'ci.low'], qtls['Fatq1', 'ci.high']), y=c(-.05,-.05), type='l',
      lwd=3, col=cols[1], pch='', lty=1, cex=1.7); points(x=qtls['Fatq1', "position.mb"],
                                           y=-.05, pch=20, col=cols[1], cex=2.7)
text(x=qtls['Fatq1', 'ci.low'], y=-.05, labels="Fatq1", pos = 2, cex=1.8, font = 3)

# Fatq2
lines(x=c(qtls['Fatq2', 'ci.low'], qtls['Fatq2', 'ci.high']), y=c(-.05,-.05), type='l',
      lwd=3, col=cols[1], pch='', lty=1, cex=1.7); points(x=qtls['Fatq2', "position.mb"],
                                           y=-.05, pch=20, col=cols[1], cex=2.7)
text(x=qtls['Fatq2', 'ci.low'], y=-.05, labels="Fatq2", pos = 2, cex=1.8, font = 3)

dev.off()



# Figure 3.pdf
# 9/26/2008 11:53:18 AM

# loads conserved tracks, and genes and coordinates
load("../data/cns.rda")
(g2<- read.delim("../data/Fatq2b_genes.txt"))
# other colors
cols = brewer.pal(5, "Set1")
colr = brewer.pal(9, "YlOrRd")
colg = brewer.pal(9, "YlGn")


pdf("../figures/Figure3Test.pdf", width = 12, height = 8)

nf<- layout(mat=matrix(c(1, 2, 2, 3, 3, 4), nrow=3, ncol=2), widths=c(1.7, 2.5), heights=c(2, 2, .7))

# plots the interval mapping and composite interval mapping for total fat mass i.e. Figure 2
par(mar=c(5,4.5,4,1)+.1, cex.lab=1.8, cex.axis = 1.6, bty='l', xaxt='s', yaxt='s', lwd = 2)    ##

plot(out.TFe.ss, col=cols[1], xlim=c(75,180), xlab='Map Position (Mb)', ylab='LOD',
     lwd=5, ylim=c(0,11))
plot(cim.TFss.w2.arg, col=cols[2], lwd=3, add=TRUE)
plot(rcim.TFss.w2.arg.400x, col='black', lwd=3, add=TRUE)
abline(h=3.44, col='gray', lwd=.7)
text(x = 80, y = 3.44, label = expression(alpha~"="~0.01), pos = 3, cex = 1.5)
mtext("A", side = 3, line = 1, outer = FALSE, at = 70, font = 2, cex = 1.8)

# blank plot for the figure legend
par(mar=c(0,4,3,0)+.1, bty='n', xaxt='n', yaxt='n', lwd = 2)

plot(rcim.TFss.w2.arg.400x, col='white', xlim=c(172,176), ylim=c(-1, 5),
     incl.markers=FALSE, show.marker.names=FALSE, bty='n',
     ylab="", xlab="", lwd=3)

legend('topleft', legend=c("Interval Mapping", "Composite Interval Mapping",
                    "CIM -400x -step=.1, -w=2 ", "",
                    "CIM -400x -step=1,  -w=0.50",
                    "CIM -400x -step=.5, -w=1.00", "",
                    "All Conserved Sequences MMU, HSA, CFA",
                    "Intergenic Conserved Sequences",
                    "Conserved Gene Elements (UTR, Exons)",
                    "miRNA"),
       lty=c(1, 1, 1, 1, 2, 3, 1, 1, 1, 1), lwd=5,
       col=c(cols[1:2], 'black', 'white', colg[8], cols[4], "white", colr[9], cols[5], "dodgerblue", cols[3]),# "chartreuse4"), 
       bty='n',
       cex=1.7)


# plots representative results of the replicated composite interval mapping
par(mar=c(5, 5, 4, 1)+.1, cex.lab=1.8, cex.axis=1.6, bty='l', xaxt='s', yaxt='s', lwd = 2)

plot(rcim.TFss.w2.arg.400x, col='black',  xlim=c(173.2, 175.7), ylim=c(-5, 7.9),
     incl.markers=TRUE, cex=1.5, ylab="Median LOD of 400 CIM Replications", xlab="Map Position (Mb)",
     show.marker.names=TRUE, lwd=4)
plot(TF.s05, lodcolum=c(3), col=cols[4], lwd=4, lty=c(3), add=TRUE)
plot(TF.s1, lodcolum=c(2), col=colg[8],lwd=4, lty=c(2), add=TRUE)
plot(cim.TFss.w2.arg, col=cols[2], lwd = 4, add = TRUE)

abline(h=3.44, col='gray', lwd=.7)
text(x = 175.59, y = 3.44, labels = expression(alpha~"="~0.01), pos = 3, cex = 2)

# plot genes and gene names under the axis
mg<- 1:length(g2[,1])
lapply(mg[-c(31,32)], function(i) {
  if (g2$strand[i]=="F"){
    arrows(x0=g2$start[i], x1=g2$end[i], y1=-g2$y[i], y0=-g2$y[i],
           angle=15, length=.12, lwd=5, col=colr[6])
    text(x=g2$end[i]+.02, y=-g2$y[i], labels=g2$symbol[i], cex=1.4, pos = 4)
  }
  
  if (g2$strand[i]=="R"){
    arrows(x0=g2$end[i], x1=g2$start[i], y1=-g2$y[i], y0=-g2$y[i],
           angle=15, length=.12, lwd=5, col=colr[6])
    text(x=g2$end[i]+ .02, y=-g2$y[i], labels=g2$symbol[i], cex=1.4, pos = 4)
    
  }
}
       
       )
axis(1, at=c(173,173.5,174,174.5, 175), labels=TRUE, tick=TRUE)
axis(2, at=c(-3, 2, 0, 1, 2, 3, 4), labels=FALSE, tick=FALSE)

mtext("B", side = 3, line = 1, outer = FALSE, at = 173.05, font = 2, cex = 1.8)

# plot sequence analysis of Fatq2b confidence interval

par(mar=c(4.2, 5, 0, 1)+.1, cex.lab=1.6, cex.axis=1.5, bty='n', xaxt='s', yaxt='n')##
plot(TF.s1, lodcolum=c(2), col = "white",  xlim=c(173.2, 175.7), ylim = c(.1, 2),
     cex=1.5, ylab="", xlab="Map Position (Mb)")
axis(1, at=c(173, 173.5, 174, 174.5, 175, 175.5, 176), labels=TRUE, tick=TRUE)

# define intervals for the conserved sequences
c3<- rbind(ce2[ce2$start > 173.440367 & ce2$end < 173.486543, ],
           ce2[ce2$start > 173.721333 & ce2$end < 173.861889, ],
           ce2[ce2$start > 173.998611 & ce2$end < 174.077834, ],
           ce2[ce2$start > 174.077747 & ce2$end < 174.110426, ],
           ce2[ce2$start > 174.267580 & ce2$end < 174.272759, ])
c3<- rbind(c3, ce2[ce2$type == "intron",], ce2[ce2$type == "exon",], ce2[ce2$type == "UTR",])

# color code 
c3$col[c3$type == "intergenic"] = cols[5]
c3$col[c3$type == "UTR"] = "dodgerblue"
c3$col[c3$type == "exon"] = "dodgerblue"
c3$col[c3$type == "intron"] = "dodgerblue"

# y axis coordinate
c3$y[c3$type == "intergenic"] = .1
c3$y[c3$type == "UTR"] = .4
c3$y[c3$type == "exon"] = .4
c3$y[c3$type == "intron"] = .4

# plotting
apply(ce2[ce2$type == "intergenic",], 1, function(i) {
  x1 = i["start"]
  x2 = i["end"]
  y1 = .20 + as.numeric(i["y"])
  y2 = y1
  lines(x = c(x1, x2), y = c(y1, y2), lwd = 6.3, col = colr[9])
})

apply(c3, 1, function(i) {
  x1 = i["start"]
  x2 = i["end"]
  y1 = .7 + as.numeric(i["y"])
  y2 = y1
  col1 = i["col"]
  lines(x = c(x1, x2), y = c(y1, y2), lwd = 6.3, col = col1)
})

apply(mirnas, 1, function(i) {
  x1 = i["start"]
  x2 = i["end"]
  y1 = 1.4 + as.numeric(i["y"])
  y2 = y1
  col1 = cols[3]
  lines(x = c(x1, x2), y = c(y1, y2), lwd = 6.3, col = col1)
})

mtext("C", side = 2, line = 0, cex = 1.8, font = 2, las = 2)

dev.off()


# Figure 4

