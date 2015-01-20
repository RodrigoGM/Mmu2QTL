## quick and dirty plots
setwd("../analysis/")
library(qtl)
source("../src/qtl.R")
source("../src/qtl2.R")

# loads output from fat.qtl.R
load("scanone_output.rda")
load("cim_output.rda")

# Interval Mapping QTL plots
# Total Fat
plot(out.TF.nc, out.TF.ss, out.TF.ss2, col=c("red", "blue","purple"), xlim=c(105,180), lty=1, lwd=2, incl.markers=TRUE, show.marker.names=FALSE, main= "QTL map for Total Fat \n Haley-Knott Scan")
legend("topright", c("No Covariates", "Sex & sac as Additive", "Sex & sac as Add and Int"), col=c('red', 'blue', 'purple'), lty=c(1,1,1), lwd=3, bty='n')                     

plot(out.TF.nc, out.TF.sex, out.TF.ss2, col=c("red", "blue", "brown"), xlim=range(105:180), lty=1, lwd=2, incl.markers=TRUE, show.marker.names=FALSE, main= "QTL map for Total Fat \n Haley-Knott Scan")
plot(out.TF.sac, out.TF.ss, col=c('purple', 'green'), xlim=range(105:180), lty=1, lwd=2, add=TRUE)
locator(n=10, type="n") 
legend('topleft', c('No Covariates', 'Sex as Additive', 'sac as Additive', 'Sex & sac as Additive', 'Sex & sac as Interactive'), col= c("red", "blue", "purple", 'green','brown'), lty=c(1,1,1,1), lwd=c(2,2,2,2), bty='n')

plot(out.TFe.nc, out.TFe.ss, out.TFe.ss2, col=c("red", "blue","purple"), xlim=range(105:180), lty=1, lwd=2, incl.markers=TRUE, show.marker.names=FALSE, main= "QTL map for Total Fat \n corrected by Strain \n Haley-Knott Method")
legend("topleft", c("No Covariates", "Sex & sac as Additive", "Sex & sac as Add and Int"), col=c('red', 'blue', 'purple'), lty=c(1,1,1), lwd=3, bty='n')

plot(out.TFe.sac, out.TFe.sex, out.TFe.ss2, col=c("purple", "blue", "brown"), xlim=range(105:180), lty=1, lwd=2, incl.markers=TRUE, show.marker.names=TRUE, main= "QTL map for Total Fat \n corrected by Strain \n Haley-Knott Scan")
plot(out.TFe.nc, out.TFe.ss, col=c('red', 'green'), xlim=range(105:180), lty=1, lwd=2, add=TRUE)
legend('topleft', c('No Covariates', 'Sex as Additive', 'sac as Additive', 'Sex & sac as Additive', 'Sex & sac as Interactive'), col= c("red", "blue", "purple", 'green','brown'), lty=c(1,1,1,1), lwd=c(2,2,2,2), bty='n')

# Individual Fat pads
plot(out.GFP.nc, out.RFP.nc, out.FFP.nc, col=c('red', 'blue', 'green'), xlim=range(143:180), lwd=2, incl.markers=TRUE, show.marker.names=TRUE, main="QTL for individual fat pads\n No Covariates")
plot(out.MFP.nc, col=c('purple', ), lty=c(1,2), add=TRUE)
locator(n=10, "l")
legend('topright', c('GFP', 'RFP', 'MFP', 'FFP'), col=c('red', 'blue', 'purple', 'green'), lty=c(1,1,1,1,2), lwd=2, bty='n')

plot(out.GFP.ss, out.RFP.ss, out.FFP.ss, col=c('red', 'blue', 'green'), xlim=range(143:180), lwd=2, incl.markers=TRUE, show.marker.names=TRUE, main="QTL for individual fat pads\n Sex and sac as Additive Covariates")
plot(out.MFP.ss, col=c('purple'), lty=c(1,2), lwd=2, add=TRUE)
locator(n=10, "l")
legend('topright', c('GFP', 'RFP', 'MFP', 'FFP'), col=c('red', 'blue', 'purple', 'green'), lty=c(1,1,1,1,2), lwd=2, bty='n')

plot(out.GFP.ss2, out.RFP.ss2, out.FFP.ss2, col=c('red', 'blue', 'green'), xlim=range(143:180), lwd=2, incl.markers=TRUE, show.marker.names=TRUE, main="QTL for individual fat pads\n Sex and sac as Additive & Interactive Covariates")
plot(out.MFP.ss2, col=c('purple'), lty=c(1,2), lwd=2, add=TRUE)
plot(out.TF.ss2, col=c('orange3'), lwd=0.5, lty=1, add=TRUE)
locator(n=10, "l")
legend('topright', c('GFP', 'RFP', 'MFP', 'FFP', "TF"), col=c('red', 'blue', 'purple', 'green', 'orange3'), lty=c(1,1,1,1,2), lwd=2, bty='n')


plot(out.GFPe.nc, out.RFPe.nc, out.FFPe.nc, col=c('red', 'blue', 'green'), xlim=range(143:180), lwd=2, incl.markers=TRUE, show.marker.names=TRUE, main="QTL for individual fat pads\n corrected by strain \n No Covariates")
plot(out.MFPe.nc, out.WK9e.nc, col=c('purple', 'black'), lty=c(1,2), add=TRUE)
locator(n=10, "l")
legend('topleft', c('GFP', 'RFP', 'MFP', 'FFP', 'WK9'), col=c('red', 'blue', 'purple', 'green', 'black'), lty=c(1,1,1,1), lwd=2, bty='n')

plot(out.GFPe.ss, out.RFPe.ss, out.FFPe.ss, col=c('red', 'blue', 'green'), xlim=range(143:180), lwd=2, incl.markers=TRUE, show.marker.names=TRUE, main="QTL for individual fat pads\n corrected by strain \n Sex and sac as Additive Covariates")
plot(out.MFPe.ss, out.WK9e.sw, col=c('purple', 'black'), lty=c(1,2), add=TRUE)
locator(n=10, "l")
legend('topleft', c('GFP', 'RFP', 'MFP', 'FFP', 'WK9'), col=c('red', 'blue', 'purple', 'green', 'black'), lty=c(1,1,1,1), lwd=2, bty='n')

plot(out.GFPe.ss2, out.RFPe.ss2, out.FFPe.ss2, col=c('red', 'blue', 'green'), xlim=range(143:180), lwd=2, incl.markers=TRUE, show.marker.names=TRUE, main="QTL for individual fat pads\n corrected by strain \n Sex and sac as Additive & Interactive Covariates")
plot(out.MFPe.ss2, out.WK9e.sw2, col=c('purple', 'black'), lty=c(1,2), add=TRUE)
plot(out.WK9.sw2, out.TF.ss2, col=c('black', 'orange'), lwd=0.5, lty=1, add=TRUE)
locator(n=10, "l")
legend('topleft', c('GFP', 'RFP', 'MFP', 'FFP', 'WK9'), col=c('red', 'blue', 'purple', 'green', 'black'), lty=c(1,1,1,1), lwd=2, bty='n')


# All Total Fat and Fat Pads
plot.mscanone(out.FATtrt.nc, title="QTL for Fat traits corrected by strain\n")
plot.mscanone(out.FATtrt.ss, title="QTL for Fat traits corrected by strain \n sex and sac wt as add", xlim=range(105:182))
plot.mscanone(out.FATtrt.ss2, title="QTL for Fat traits corrected by strain \n sex and sac wt as add and int\n")

plot.mscanone(out.FATtrt.ss, title="QTL for Fat traits \n corrected by strain, sex and sac wt\n", xlim=range(105:182), lwd=3)
abline(h=2.27, col="black", lty=2, lwd=2)
text(x=105, y=2.8, labels="0.05", cex=1.1)


# qtl and effect plots
nf <- layout(matrix(c(1,2,1,3), 2, 2, byrow=TRUE), c(3,1), respect=TRUE)
layout.show(nf)
par(mar=c(0,4,0,.5))
plot.mscanone(out.FATtrte.ss, xlim=range(105:180), lwd=4, ylab="LOD", xlab="Map Position (Mbp)", col=c('red', 'blue', 'brown', 'green', 'purple'), mtitle=FALSE, title="", incl.legend=FALSE, incl.markers=FALSE)
abline(h=2.27, col="black", lty=2, lwd=3)
#abline(h=0, col="black", lty=1, lwd=1)
text(x=105, y=2.8, labels="0.05", cex=1.5)
text(x=163, y=13, labels="Fatq2", cex=1.5)
text(x=105, y=13.5, labels="[A]", cex=2)
legend(x=102,y=13,  c('TF', 'GFP', 'RFP', 'MFP', 'FFP'), col=c('red', 'blue', 'brown', 'green', 'purple'), lty=c(1,1,1,1), lwd=4.25, bty='n', cex=1)
par(mar=c(1,4,0,1))
rs27669516<- find.marker(D2sq, chr=2, pos=174)
effectplot(D2sq, pheno.col=find.pheno(D2sq, "TF.ss"), mname1=rs27669516, main="", ylab="Total Fat (g)")
text(x=3, y=0.021, labels="[B]", cex=2)
text(x=1, y=-.08, labels="174.0 Mb", cex=1, adj=c(.1,0))
par(mar=c(0,4,2,1))
rs27307998<- find.marker(D2sq, chr=2, pos=156.97617)
effectplot(D2sq, pheno.col=find.pheno(D2sq, "TF.ss"), mname1=rs27307998, main="", ylab="Total Fat (g)")
text(x=3, y=0.02, labels="[C]", cex=2)
text(x=1, y=-.062, labels="156.9 Mb", cex=1, adj=c(.1,0))

# HG2D plots
plot.mscanone(Fat.ss.2D, title="QTL for Fat traits \n HG2D uncorrected data, sex and sac", xlim=range(139:181))
plot.mscanone(Fate.ss.SBC, title="QTL for Fat traits \n SBC data corrected by strain, sex and sac wt", xlim=range(139:181))
plot.mscanone(out.FATtrt.ss, title="QTL for Fat traits \n corrected by strain, sex and sac wt", xlim=range(139:181))

plot.mscanone(Fat.ss.2D.em, title="QTL for Fat tratis\nHG2D uncorrected data, sex and sac, EM", xlim=range(139:181))


# CIM plots
# plotting interval mapping, cim, and rcim
plot(out.TFe.ss, cim.TFss.w2.arg, rcim.TFss.w2.imp.400x, xlim=range(168:179), lwd=2, ylab="LOD", col=c('red', 'blue', 'black'), incl.markers=TRUE, show.marker.names=TRUE)
abline(h=2.27, col="black", lty=3, lwd=1.5)
text(x=168.5, y=2.8, labels="0.05", cex=1)
text(x=163, y=10, labels="Fatq2", cex=1.5)
#text(x=176, y=10, labels="[B]", cex=2, adj=c(.75,.5))
legend("topright",  c('TF: Interval Mapping', 'TF: CIM 1x', 'TF: CIM 400x'), col=c('red', 'blue', 'black'), lty=1, lwd=2.2, bty='n', cex=1)
abline(v=c(170.3, 173.25, 175.5), col=c("gray", "black","black", lty=3, lwd=1))

plot(TFss.rcim, xlim=range(174:175), main="Total Fat composite interval mapping", ylab="LOD", incl.markers=TRUE, show.marker.names=TRUE)

plot(cim.TFss.w2.arg, rcim.TFss.w2.imp.50x, rcim.TFss.w2.arg.50x, col=c("red", "blue","purple"), xlim=range(105:180), lty=1, lwd=2, incl.markers=TRUE, show.marker.names=FALSE, main= "QTL map for Total Fat; HK CIM Scan")
legend("topleft", c("CIM W=2 S=.2 1x", "CIM W=2 S=.2 50x imp", "CIM W=2 S=.2 50x arg"), col=c('red', 'blue', 'purple'), lty=c(1,1,1), lwd=3, bty='n')                     

plot(cim.TFss.w2.arg, rcim.TFss.w2.imp.50x, rcim.TFss.w2.arg.50x, col=c("red", "blue","purple"), xlim=range(171:176.5), lty=1, lwd=2, incl.markers=TRUE, show.marker.names=FALSE, main= "QTL map for Total Fat; HK CIM Scan")
legend("topleft", c("CIM W=2 S=.2 1x", "CIM W=2 S=.2 50x imp", "CIM W=2 S=.2 50x arg"), col=c('red', 'blue', 'purple'), lty=c(1,1,1), lwd=3, bty='n')                     

plot(cim.TFss.w2.arg, rcim.TFss.w2.imp.50x, rcim.TFss.w2.arg.50x, col=c("red", "blue","purple"), xlim=range(148:152), lty=1, lwd=2, incl.markers=TRUE, show.marker.names=FALSE, main= "QTL map for Total Fat; HK CIM Scan")
legend("topleft", c("CIM W=2 S=.2 1x", "CIM W=2 S=.2 50x imp", "CIM W=2 S=.2 50x arg"), col=c('red', 'blue', 'purple'), lty=c(1,1,1), lwd=3, bty='n')                     

## pdf(file="../Marker CIM lod histograms.pdf")
sapply(markers, function(m) {hist(as.numeric(rcim.results[m,3:52]), xlab=m, main="Total Fat LOD scores m=3, w=1 50x CIM ", breaks=200)})
## dev.off()


nf <- layout(matrix(c(1,2,1,0), 2, 2, byrow=TRUE), c(4,3), c(2,2), respect=TRUE)
layout.show(nf)
##PLOT [A]
par(mar=c(4.1,4,0,1))
plot(out.TFe.ss, cim.TFss.w2.arg, rcim.TFss.w2.imp.400x, xlim=range(105:180), lwd=2, ylab="LOD", col=c('red', 'blue', 'black'))
abline(h=2.27, col="black", lty=3, lwd=1.5)
text(x=106.8, y=2.8, labels="0.05", cex=1.1)
text(x=163, y=10, labels="Fatq2", cex=1.5)
text(x=106, y=10, labels="[A]", cex=2, adj=c(0,.50))
legend(x=102,y=9,  c('TF: Interval Mapping', 'TF: CIM 1x', 'TF: CIM 400x'), col=c('red', 'blue', 'black'), lty=1, lwd=3, bty='n', cex=1.2)
##PLOT [B]
par(mar=c(0,1,0,1))
plot(out.TFe.ss, cim.TFss.w2.arg, rcim.TFss.w2.imp.400x, xlim=range(168:179), lwd=2, ylab="LOD", col=c('red', 'blue', 'black'), incl.markers=TRUE, show.marker.names=TRUE)
abline(h=2.27, col="black", lty=3, lwd=1.5)
text(x=168.5, y=2.8, labels="0.05", cex=1)
#text(x=163, y=10, labels="Fatq2", cex=1.5)
text(x=176, y=10, labels="[B]", cex=2, adj=c(.75,.5))
#legend(x=173,y=9,  c('TF: Interval Mapping', 'TF: CIM 1x', 'TF: CIM 400x'), col=c('red', 'blue', 'black'), lty=1, lwd=2.2, bty='n', cex=1)

