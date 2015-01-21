###########################################################
##  These functions were writen by Ricardo Verdugo Ph.D  ##
##  at the University of California at Davis             ##
##  the functions can be freely "AS IS" without any      ##
##  warranty.                                            ##
###########################################################
countrecs=function(genotypes, missing=NA){
  if(!is.null(dim(genotypes))) { stop("Argument 'genotypes' should be a vector") }
  if(!(length(genotypes)>1)){  stop("Vector 'genotypes' hast to have more than one value")}
  count=0
  prev=NULL
  if(is.numeric(genotypes)) {
    for (i in 1:length(genotypes)){
        if(!is.na(missing) & genotypes[i]==missing){
                genotypes[i]=NA
        }  
        if(is.na(genotypes[i])){
            next
        }
        if(is.null(prev)){
          prev=genotypes[i]
          next
        } else {
          if(prev!=genotypes[i]){
            count=count+abs(genotypes[i]-prev)
          }
        }
    prev=genotypes[i]
    }
  } else {
    for (i in 1:length(genotypes)){
        if(!is.na(missing) & genotypes[i]==missing){
                genotypes[i]=NA
        }  
        if(is.na(genotypes[i])){
            next
        }
        if(is.null(prev)){
          prev=genotypes[i]
          next
        } else {
          if(prev!=genotypes[i]){
            count=count+1
          }
        }
    prev=genotypes[i]
   }
  }
    return(count)
}  # countrecs               

# Plot Makers Effect On Two contiguous intervals
# Description
# Plot function for effect sizes of three markers in two contiguous intervals
# Usage
# plot.inttest(cross, pheno.name, chr, marker, title=NULL, ...)
#
# Arguments
#       cross        a cross object
#       pheno.name   character string of phenotype name
#       chr          character string of chromosome name
#       marker       character string of marker name
#       title        character string used to title the figure
#       ...          furher arguments passed to the plot function

plot.inttest=function(cross, pheno.name, chr, marker) {
  flanking=find.flanking.markers(cross, chr, marker)
  recs=apply(cross$geno[[chr]]$data[,25:28], 1, countrecs)>0
  effectplot(recs25_28.subset, pheno.col=find.pheno(cross, "MFP"), mname1="rs13482893", main="MFP:rs13482893\nin rs33763519-rs33779311 recs", ylim=c(.165, .261))

} # plot.inttest

# find.flanking.markers finds the markers flanking to a given marker
find.flanking.markers=function(cross, chr, marker) {
  "cross" %in% class(cross) || stop("Argument cross should be of class cross (see scanone).")
  require(qtl) || stop("R/QTL package not vailable.")
  markersvec=names(pull.map(cross)[[1]])
  marker %in% markersvec || stop(paste("Marker", marker, "could not be found in the map."))
  markeroffset=match(marker, markersvec)
  leftflank=character(0)
  if(markeroffset>1){
    leftflank=markersvec[markeroffset-1]
    } else {
    leftflank=NA
  }
  rightflank=character(0)
    if(markeroffset<length(markersvec)){
      rightflank=markersvec[markeroffset+1]
    } else {
      rightflank=NA
    }
  return(c(leftflank, rightflank))
} # find.flanking.markers

# find.marker.pos finds the position of a marker in a cross
find.marker.pos=function(cross, chr, marker){
  "cross" %in% class(cross) || stop("Argument cross should be of class cross (see scanone).")
  require(qtl) || stop("R/QTL package not vailable.")
  marker %in% names(pull.map(cross)[[1]]) || stop(paste("Marker", marker, "could not be found in the map."))

  position=pull.map(cross)[[1]][marker]
  return(position)
} # find.marker.pos

mident=function(scan, chr) {
  if(!class(scan)[[1]]=="scanone") {
    stop("Argument scan needs to be of class 'scanone'")
  }
  identify(scan[scan[,1]==chr, 2:3], labels=rownames(scan[scan[,1]==chr,]))
 }
 

inttest=function(cross, chr, pheno.name, marker, n=3, alpha=0.4, add=F) {
# Implements a step-wise test to find the smallest
# marker interval harboring a QTL
  require("qtl") || stop("qtl package not available.")
    if(n<1) stop("A positive integer must be input for argument 'n'.")
  map=pull.map(cross)[[chr]]
  genotable=cross$geno[[chr]]$data
  phenotypes=colnames(cross$pheno)
  # Store phenotypes
    pheno=cross$pheno[,match(pheno.name, phenotypes)]

  if(!pheno.name %in% phenotypes ) { stop(paste(pheno.name, "not found in cross.")) }
  if(is.na(match(marker, names(map)))){
    stop(paste("Marker", marker, "couldn\'t be found in cross"))
  }
  M2=match(marker, names(pull.map(cross)[[chr]]))
  if(M2==1)M2=2

  int.obj=create_int(pheno, genotable, M2, n)
  # Create subset of recombinant individuals
  #recs=apply(genotable[,M1:M2], 1 ,countrecs)>0

  print("## First step ##")
  # Test markers
  M_pval=M_test(int.obj, add)

  if(M_pval[2]<alpha) {
    print("P2: yes")
    M2=M2+1 # move o the right
  } else {
    print("P2: no")
    if(M_pval[1]<alpha) {
      print("P1: yes")  
      M2=attr(int.obj,"M1") # move to the left
    } else {
      print("P2: no")
        M2=M2+1 # move o the right
    }
  }
  print(paste("M1=", attr(int.obj, "M1")))
  print(paste("M2=", attr(int.obj, "M2")))
  
  ## Following steps
  out_left=NULL
  out_right=NULL
  int.p.obj=int.obj
  nstep=0
  
  while (is.null(out_left) | is.null(out_right)) {
    nstep=nstep+1
    print(paste("## Step", nstep, "##"))
    # Update interval object
    int.obj=create_int(pheno, genotable, M2, n)
    print(paste("M1=", attr(int.obj, "M1")))
    print(paste("M2=", attr(int.obj, "M2")))
    # Test markers
    M_pval=M_test(int.obj, add)
    if(M_pval[2]<alpha) {
      print("P2: yes")    
      M2=M2+1 # move to the right
    } else {
      print("P2: no")
      if(M_pval[1]<alpha) {
        print("P1: yes")      
        out_left=attr(int.p.obj,"M1") # return left limit from previous run. Q2 is true
        out_right=M2
      } else { # Q3 is true
        print("P1: no")
        out_left=attr(int.obj,"M1")
        out_right=M2
      }                                          
    }
  } # while
  lim_markers=map[c(out_left, out_right)]
  return(lim_markers)
} # inttest

M_test=function(int_obj, add=F) {
  if(!add){
    indep1=as.factor(int_obj[,2])
    indep2=as.factor(int_obj[,3])
  } else {
    indep1=int_obj[,2]
    indep2=int_obj[,3]
  }

  lm.obj1=anova(lm(pheno~indep1, data=int_obj))$"Pr(>F)"
  lm.obj2=anova(lm(pheno~indep2, data=int_obj))$"Pr(>F)"
  out=c(lm.obj1[1], lm.obj2[1])
  names(out)=c("M1","M2")
  print(out)
  return(out)

} # M_test

M_test_model=function(int_obj, add=F) {
  n=nrow(int_obj)
  x1.1=int_obj[,2]
  x2.1=int_obj[,3]
  x1.2=int_obj[,2]==median(unique(int_obj[,2]), na.rm=TRUE)
  x2.2=int_obj[,3]==median(unique(int_obj[,3]), na.rm=TRUE)
  x1.formula=NULL
  x2.formula=NULL
  outrownames=character(0)
  
  n.factors=0
  if(!add){
    x1.formula=pheno~x1.1 + x1.2
    x2.formula=pheno~x2.1 + x2.2
    n.factors=3
    outrownames=c("alpha", "delta", "pval{a}", "pval{d}")
  } else {
    x1.formula=pheno~x1.1
    x2.formula=pheno~x2.1
    n.factors=2
    outrownames=c("alpha", "pval{a}")

  }
  codeddata=data.frame(pheno=int_obj$pheno, x1.1, x1.2, x2.1, x2.2) 
  lm.obj1=lm(x1.formula, data=codeddata)
  lm.obj2=lm(x2.formula, data=codeddata)
  anova.obj1=summary(lm.obj1)
  anova.obj2=summary(lm.obj2)
 # print(lm.obj1)
 # print(lm.obj2)
  out=cbind(lm.obj1$"coefficients"[2:n.factors], lm.obj2$"coefficients"[2:n.factors])
  out=rbind(out, cbind(anova.obj1$coefficients[2:n.factors,"Pr(>|t|)"], anova.obj2$coefficients[2:n.factors,"Pr(>|t|)"]))
  colnames(out)=c("M1","M2")
  rownames(out)=outrownames
  attr(out, "M1")=attr(int_obj, "M1")
  attr(out, "M2")=attr(int_obj, "M2")
  attr(out, "Min.n")=attr(int_obj, "Min.n")
  attr(out, "M1.name")=attr(int_obj, "M1.name")
  attr(out, "M2.name")=attr(int_obj, "M2.name")
  return(out)
} # M_test


old_M_test=function(pheno, recs, genotable, M1, add) {
   indep=genotable[recs,M1]
   if(!add) indep=as.factor(indep)
   lm.obj=anova(lm(pheno[recs]~indep))$"Pr(>F)"
   return(lm.obj[1])
} # left_test

minnint=function(cross, chr, pheno, genotable, M2, n) {
  # This function finds the closests falnking markers to 
  # a middle marker allowing a minimum number recombinant
  # individuals per genotype of flanking marker.
  
  require("qtl") || stop("qtl package not available.")
  M1=M2
  celln=0
  while (celln<n) {
    # Check if end of chromosome
    if(M1==0) stop("Left end of chromosome reached. Posible insufficient recombinants")
    #print(paste(celln, "vs", n))
    # Update left marker
    M1=M1-1
    # Find recombinant individuals
    recs=apply(genotable[,M1:M2], 1 ,countrecs)>0
    # Create subset
    subset=subset(cross, ind=recs)
    # Create factor variables
    fact1=as.factor(genotable[, M1])
    #fact2=as.factor(genotable[,M2])
    # Count number of non NA values per left
    # marker genotype in rec subset
    celln=min(table(fact1[recs & !is.na(pheno)]))
    #print(table(fact1[recs & !is.na(pheno)]))
    #print(paste("End loop:,"\n","celln:", celln))
    #print(paste("M1:", M1))
  }
  return(M1)
} # minnint

create_int=function(pheno, genotable, M2, n) {
  # This function finds the closests falnking markers to 
  # a middle marker allowing a minimum number recombinant
  # individuals per genotype of flanking marker.
  
  require("qtl") || stop("qtl package not available.")
  nmin=5 # in number of recombinnat animals
  M1=M2
  celln=0
  recs=logical(0)
  while (celln<n) {
    # Check if end of chromosome
    if(M1==0) stop("Left end of chromosome reached. Posible insufficient recombinants")
    #print(paste(celln, "vs", n))
    # Update left marker
    M1=M1-1
    # Find recombinant individuals
    recs=apply(genotable[,M1:M2], 1 ,countrecs)>0
    if(sum(recs)<nmin) next
    # Create factor variables
    fact1=as.factor(genotable[, M1])
    #fact2=as.factor(genotable[,M2])
    # Count number of non NA values per left
    # marker genotype in rec subset
    celln=min(table(fact1[recs & !is.na(pheno)]))
    #print(table(fact1[recs & !is.na(pheno)]))
    #print(paste("End loop:,"\n","celln:", celln))
    #print(paste("M1:", M1))
  }
  pheno=pheno[recs]
  M1_geno=genotable[recs, M1]
  M2_geno=genotable[recs, M2]
  output=data.frame(pheno, M1_geno, M2_geno)
  attr(output, "M1")=M1
  attr(output, "M2")=M2
  attr(output, "Min.n")=celln
  return(output)
} # create_int

old_create_int=function(cross, chr, pheno, genotable, M2, n) {
  # This function finds the closests falnking markers to 
  # a middle marker allowing a minimum number recombinant
  # individuals per genotype of flanking marker.
  
  require("qtl") || stop("qtl package not available.")
  nmin=5 # in number of recombinnat animals
  M1=M2
  celln=0
  recs=logical(0)
  while (celln<n) {
    # Check if end of chromosome
    if(M1==0) stop("Left end of chromosome reached. Posible insufficient recombinants")
    #print(paste(celln, "vs", n))
    # Update left marker
    M1=M1-1
    # Find recombinant individuals
    recs=apply(genotable[,M1:M2], 1 ,countrecs)>0
    if(sum(recs)<nmin) next
    # Create subset
    subset=subset(cross, ind=recs)
    # Create factor variables
    fact1=as.factor(genotable[, M1])
    #fact2=as.factor(genotable[,M2])
    # Count number of non NA values per left
    # marker genotype in rec subset
    celln=min(table(fact1[recs & !is.na(pheno)]))
    #print(table(fact1[recs & !is.na(pheno)]))
    #print(paste("End loop:,"\n","celln:", celln))
    #print(paste("M1:", M1))
  }
  pheno=pheno[recs]
  M1_geno=genotable[recs, M1]
  M2_geno=genotable[recs, M2]
  output=data.frame(pheno, M1_geno, M2_geno)
  attr(output, "M1")=M1
  attr(output, "M2")=M2
  attr(output, "Min.n")=celln
  return(output)
} # create_int

create_side_int=function(pheno, genotable, M2, n, side="left") {
  # This function finds the closests falnking markers to 
  # a middle marker allowing a minimum number recombinant
  # individuals per genotype of flanking marker.
  
  require("qtl") || stop("qtl package not available.")
  nmin=5 # in number of recombinnat animals
  M1=M2
  celln=0
  recs=logical(0)
  while (celln<n) {
    # Check if end of chromosome
    if(M1==0) stop("Left end of chromosome reached. Posible insufficient recombinants")
    #print(paste(celln, "vs", n))
    # Movement depends on side
    if(side=="right") {
      M2=M2+1 # move to the right
    } else {
      M1=M1-1 #move to the left
    }
    # Find recombinant individuals
    recs=apply(genotable[,M1:M2], 1 ,countrecs)>0
    if(sum(recs)<nmin) next
    # Create factor variable
 #   if(side=="right") {
       fact1=as.factor(genotable[, M1])
  #  } else {
       fact2=as.factor(genotable[, M2])
  #  }
    # Count number of non NA values per left
    # marker genotype in rec subset
    celln=min(table(fact1[recs & !is.na(pheno)]), table(fact2[recs & !is.na(pheno)]))
    #print(table(fact1[recs & !is.na(pheno)]))
    #print(paste("End loop:,"\n","celln:", celln))
    #print(paste("M1:", M1))
  }
  pheno=pheno[recs]
  M1_geno=genotable[recs, M1]
  M1.name=colnames(genotable)[M1]
  M2_geno=genotable[recs, M2]
  M2.name=colnames(genotable)[M2]
  output=data.frame(pheno, M1_geno, M2_geno)
  attr(output, "M1")=M1
  attr(output, "M2")=M2
  attr(output, "Min.n")=celln
  attr(output, "M1.name")=M1.name
  attr(output, "M2.name")=M2.name
  return(output)
} # create_side_int

cross2int_3pt=function(cross, chr, pheno, marker, n) {
  require("qtl") || stop("qtl package not available.")
  phenoN=find.pheno(cross, pheno)
  phenovec=cross$pheno[,phenoN]
  genotab=cross$geno[[chr]]$data
  if(is.numeric(marker))
    markerN=marker
  else
    markerN=match(marker, colnames(genotab))

  if(length(markerN)>1){
    int_3pt=lapply(markerN, function(x) create_int_3pt(phenovec, genotab, x, n))
    for (i in 1:length(int_3pt)) attr(int_3pt[[i]], "phenotype")=pheno
  } else {
    int_3pt=create_int_3pt(phenovec, genotab, markerN, n)
    attr(int_3pt, "phenotype")=pheno
    }
  return(int_3pt)
} # cross2int_3pt


create_int_3pt=function(pheno, genotable, M2, n) {

  left_side=create_side_int(pheno, genotable, M2, n, "left")
  right_side=create_side_int(pheno, genotable, M2, n, "right")
  output=list(left=left_side, right=right_side)
  class(output)=c("list", "int_3pt")
  return(output)
} # create_3pt_int

plot.int_3pt=function(int_3pt, add=F, ...) {
  par(mfrow=c(2,2), mar=c(6.1, 4.1, 2, 1))
  means=geno.means.int_3pt(int_3pt)  
  ses=geno.se.int_3pt(int_3pt)
  phenoname=attr(int_3pt, "phenotype")
  ylims=range(c(as.vector(means+ses), as.vector(means-ses)))
  xlims=c(0.5,3.5)
  tests=M_test_3pt(int_3pt, add)
  M1=attr(tests, "M1.name")
  M2=attr(tests, "M2.name")
  M3=attr(tests, "M3.name")  
  M=c(M1, M2, M2, M3)
  lw=diff(xlims)*.02
  pw=diff(ylims)*.02

  for(i in 1:4) {
    xlabel=paste(M[i], "\n", "a=", round(tests[i,"alpha"], 2), " (", format.pval(tests[i,"pval{a}"], digits=3, eps=1e-3), ")", sep="")
    if("delta" %in% colnames(tests))
      xlabel=paste(xlabel, paste("d=", round(tests[i,"delta"], 2), " (", format.pval(tests[i,"pval{d}"], digits=3, eps=1e-3), ")", sep=""), sep="\n")
    plot(1:3, means[i,], type="b", xlab="", ylab=phenoname, xlim=xlims, ylim=ylims, frame.plot=TRUE, axes=F, 
       main=paste(phenoname, ": ", M[i], sep=""))
    mtext(xlabel, 1, 2, padj=1)
    axis(1, 1:3, labels=c("AA", "AB", "BB"))
    axis(2, at=round(seq(ylims[1], ylims[2], by=diff(ylims)/4), 2), las=2)
    arrows(1:3, means[i,]+pw, 1:3, means[i,]+ses[i,], angle=90, lty=2, code=0)
    segments((1:3)-lw, means[i,]+ses[i,], (1:3)+lw, means[i,]+ses[i,])
    arrows(1:3, means[i,]-pw, 1:3, means[i,]-ses[i,], angle=90, lty=2, code=0)
    segments((1:3)-lw, means[i,]-ses[i,], (1:3)+lw, means[i,]-ses[i,])

  }
} # plot.test_3pt

geno.se.int_3pt=function(int_3pt) {
# internal function
  se=function(x) sd(x, na.rm=T)/sqrt(length(x[!is.na(x)]))
  matrix(c(tapply(int_3pt$left$pheno, int_3pt$left$M1_geno, se), 
           tapply(int_3pt$left$pheno, int_3pt$left$M2_geno, se),
           tapply(int_3pt$right$pheno, int_3pt$right$M1_geno, se),
           tapply(int_3pt$right$pheno, int_3pt$right$M2_geno, se)), nrow=4, byrow=T) 
} # geno.se.int_3pt

geno.means.int_3pt=function(int_3pt) {
# Internal function
  matrix(c(tapply(int_3pt$left$pheno, int_3pt$left$M1_geno, mean, na.rm=T), 
           tapply(int_3pt$left$pheno, int_3pt$left$M2_geno, mean, na.rm=T),
           tapply(int_3pt$right$pheno, int_3pt$right$M1_geno, mean, na.rm=T),
           tapply(int_3pt$right$pheno, int_3pt$right$M2_geno, mean, na.rm=T)), nrow=4, byrow=T)
}

M_test_3pt=function(int_3pt, add=F) {
  left_test=M_test_model(int_3pt$"left", add)
  right_test=M_test_model(int_3pt$"right", add)  
  output=cbind(left_test, right_test)
  N=N.int_3pt(int_3pt)
  n=n.int_3pt(int_3pt)
  mean=mean.int_3pt(int_3pt)
  se=se.int_3pt(int_3pt)
  output=rbind(N, n, mean, se, test=output)
  colnames(output)=c("M1", "M2", "M2.2", "M3")
  output=as.data.frame(as.matrix(t(output)))
  attr(output, "M1")=attr(int_3pt[[1]], "M1")
  attr(output, "M2")=attr(int_3pt[[1]], "M2")
  attr(output, "M2.2")=attr(int_3pt[[2]], "M1")
  attr(output, "M3")=attr(int_3pt[[2]], "M2")
  attr(output, "Left-Min.n")=attr(int_3pt[[1]], "Min.n")
  attr(output, "Right-Min.n")=attr(int_3pt[[2]], "Min.n")
  attr(output, "M1.name")=attr(int_3pt[[1]], "M1.name")
  attr(output, "M2.name")=attr(int_3pt[[1]], "M2.name")
  attr(output, "M3.name")=attr(int_3pt[[2]], "M2.name")
  class(output)=c("test_3pt", "data.frame")
  return(output)

} # M_test

mean.int_3pt=function(int_3pt){
 left_mean=apply(int_3pt[[1]][,2:3],2, mean, na.rm = T)
 right_mean=apply(int_3pt[[2]][,2:3],2, mean, na.rm = T)
 output=c(left_mean, right_mean)
 return(output)
} # mean.int_3pt

N.int_3pt=function(int_3pt) {
 left_N=apply(!is.na(int_3pt[[1]][,1]) & !is.na(int_3pt[[1]][,2:3]),2, sum)
 right_N=apply(!is.na(int_3pt[[2]][,1]) & !is.na(int_3pt[[2]][,2:3]),2, sum)
 output=c(left_N, right_N)
 return(output)
} # N.inst_3pt

n.int_3pt=function(int_3pt) {
  left_n=numeric(2)
  right_n=numeric(2)
  left_n[1]=min(table(int_3pt[[1]][!is.na(int_3pt[[1]][,1]),2]))
  left_n[2]=min(table(int_3pt[[1]][!is.na(int_3pt[[1]][,1]),3]))
  right_n[1]=min(table(int_3pt[[2]][!is.na(int_3pt[[2]][,1]),2]))
  right_n[2]=min(table(int_3pt[[2]][!is.na(int_3pt[[2]][,1]),3]))
 output=c(left_n, right_n)
 return(output)
} # mean.int_3pt

se.int_3pt=function(int_3pt) {
  se=function(x) sd(x, na.rm=T)/sqrt(length(x[!is.na(x)]))
  left_se=apply(int_3pt[[1]][,2:3],2, se)
  right_se=apply(int_3pt[[2]][,2:3],2, se)
  output=c(left_se, right_se)
  return(output)
} # se.int_3pt

old.summary.test_3pt=function(test_3pt, eps=1e-15) {
 coln=20
 colsep=5
 pvalwd=7
 rowsnames=rownames(test_3pt)
 nrows=nrow(test_3pt)
 markers=c(attr(test_3pt, "M1.name"), attr(test_3pt, "M2.name"), attr(test_3pt, "M3.name"))
 positions=c(attr(test_3pt, "M1"), attr(test_3pt, "M2"), attr(test_3pt, "M3"))
 pvals=format.pval(test_3pt[(nrows-1):nrows,], eps=eps)
 test_3pt=test_3pt[-((nrows-1):nrows),]
 pvals=matrix(pvals, ncol=4, nrow=2)
 rownames(pvals)=rowsnames[(nrows-1):nrows]
 markline=sprintf("         %10s         %10s        %1s", markers[1], markers[2], markers[3])
 posl=sprintf("Marker #     %3i              %3i               %-3i", positions[1], positions[2], positions[3])
 ruler="               +------x------+    +------x------+"
 matl1=sprintf("%-7s %8s    %10s    %-10s    %-10s", rowsnames[1], test_3pt[1,1], test_3pt[1,2], test_3pt[1,3], test_3pt[1,4])
 matl2=sprintf("%-7s %8s    %10s    %-10s    %-10s", rowsnames[2], test_3pt[2,1], test_3pt[2,2], test_3pt[2,3], test_3pt[2,4])
 matl=numeric(0)
 for(i in 3:nrow(test_3pt)) {
  matl=cbind(matl,sprintf("%-7s %8.1f    %10.1f    %-10.1f    %-10.1f", rowsnames[i], test_3pt[i,1], test_3pt[i,2], test_3pt[i,3], test_3pt[i,4]))
 }
 matl7=sprintf("%-7s %8s    %10s    %-10s    %-10s", rownames(pvals)[1], pvals[1,1], pvals[1,2], pvals[1,3], pvals[1,4])
 matl8=sprintf("%-7s %8s    %10s    %-10s    %-10s", rownames(pvals)[2], pvals[2,1], pvals[2,2], pvals[2,3], pvals[2,4]) 
 cat("", markline, posl, ruler, matl1, matl2, matl[1], matl[2], matl7, matl8, sep="\n")
 #print(outtext)
} # summary.3pt_int

summary.test_3pt=function(test_3pt, eps=1e-15) {
 D=FALSE
 if(length(test_3pt)>6) { # delta included in the model
   D=TRUE
 }
 markers=c(attr(test_3pt, "M1.name"), attr(test_3pt, "M2.name"), attr(test_3pt, "M3.name"))
 positions=c(attr(test_3pt, "M1"), attr(test_3pt, "M2"), attr(test_3pt, "M3"))
 markline=sprintf("         %10s         %10s        %1s", markers[1], markers[2], markers[3])
 posl=sprintf("Marker #     %3i              %3i               %-3i", positions[1], positions[2], positions[3])
 ruler="               +------x------+    +------x------+"
 matl=sprintf("%-7s %8s    %10s    %-10s    %-10s", "N", test_3pt$N[1], test_3pt$N[2], test_3pt$N[3], test_3pt$N[4])
 matl=c(matl, sprintf("%-7s %8s    %10s    %-10s    %-10s", "n", test_3pt$n[1], test_3pt$n[2], test_3pt$n[3], test_3pt$n[4]))
 matl=c(matl, sprintf("%-7s %8.2f    %10.2f    %-10.2f    %-10.2f", "mean", test_3pt$mean[1], test_3pt$mean[2], test_3pt$mean[3], test_3pt$mean[4]))
 matl=c(matl, sprintf("%-7s %8.2f    %10.2f    %-10.2f    %-10.2f", "se", test_3pt$se[1], test_3pt$se[2], test_3pt$se[3], test_3pt$se[4]))
 matl=c(matl, sprintf("%-7s %8.3f    %10.3f    %-10.3f    %-10.3f", "alpha", test_3pt$alpha[1], test_3pt$alpha[2], test_3pt$alpha[3], test_3pt$alpha[4]))
 if(D) {
   matl=c(matl, sprintf("%-7s %8.3f    %10.3f    %-10.3f    %-10.3f", "delta	", test_3pt$delta[1], test_3pt$delta[2], test_3pt$delta[3], test_3pt$delta[4]))
 }
 matl=c(matl, sprintf("%-7s %8.5g    %10.5g    %-10.5g    %-10.5g", "pval{a}", test_3pt$"pval{a}"[1], test_3pt$"pval{a}"[2], test_3pt$"pval{a}"[3], test_3pt$"pval{a}"[4]))
 if(D) {
   matl=c(matl, sprintf("%-7s %8.5g    %10.5g    %-10.5g    %-10.5g", "pval{d}", test_3pt$"pval{d}"[1], test_3pt$"pval{d}"[2], test_3pt$"pval{d}"[3], test_3pt$"pval{d}"[4]))
 }
 cat("", markline, posl, ruler, matl, sep="\n")
 #print(outtext)
} # summary.3pt_int

####
# Import QTL Cartographer scanning output
old.import.qtlcart=function(filename) {
  output=read.table(filename, sep="", row.names=NULL, header=T, skip=21, blank.lines.skip=F, comment.char="#", as.is=T)
  print(dim(output))
  # output=output[!output[,1]=="",]
 # print(length(grep("[[:alpha:][:punct:]]",output)))
  #output=output[-grep("[[:alpha:][:punct:]]",output)]
    print(dim(output))  
  output=as.data.frame(lapply(output, as.numeric))
 #   print(dim(output))
  return(output)
} # import.qtlcart

import.qtlcart=function(filename, method=NULL) {
  ismap=0
  output=numeric(0)
  headers=character(0)
  outattrs=list()
  map=list()
  lines=scan(filename, character(), sep="\n", quiet=TRUE)
  for(line in lines) {
    ## trim trailing white space
    line=sub('[[:space:]]+$', '', line) 
    line=sub('^[[:space:]]+', '', line)
    
    if(ismap) {
      chrmap=str2map(line)
      #map[[attr(chrmap, "Chr")]]=chrmap
      # using chr number instead of label because this si used for scaning result in cartographer.
      map[[as.character(ismap)]]=chrmap
      ismap=0
    }
    if(length(grep("^#", line))>0) {
      next
    }
    param=grep("^-", line, value=T)
    if (length(param)>0) { # This is a paramented
       if(line=="-s" | line=="-e") next # beginning or end of data
       #print(param)
       pvpair=extrpar(param)
       if(pvpair[1]=="cnum") {
         ismap=as.numeric(pvpair[2])
         next
       }
       outattrs[[pvpair[1]]]=pvpair[-1]
      #       attr(output, pvpair[1])=pvpair[-1] # this is a parameter from QTL cartographer
       next
    } # if param
    hasletters=grep("[[:alpha:]]", line, value=T)
    if (length(hasletters)>0) {
      if(length(grep("^c ", line))>0) {
        headers=unlist(strsplit(line, "\\s+", perl=T)) # headers
        next
      } else {
        next # skip lines with alphabetic characters
      }
    }
    numbers=unlist(strsplit(line, "\\s+", perl=T))
    outline=as.numeric(numbers)
    output=rbind(output, outline)
  } # for
  colnames(output)=headers
  rownames(output)=1:nrow(output)
  output=as.data.frame(output)
  output=assignattrs(output, outattrs)
  attr(output, "map")=map
  attr(output, "method")=method
  class(output)=c("scanone", "data.frame")
  return(output)
} # import.qtlcart
  
  # subroutines
  extrpar=function(str) {
    parlist=unlist(strsplit(str, "[[:space:]]+"))
    parlist[1]=sub('^-','', parlist[1])
    if(parlist[1] %in% c("window", "background", "Model", "trait", "cross")) 
      parlist=parlist[1:2] # eliminate comments from these paramenters
      #print(parlist)
    return(parlist)
  } # extrpar
    
  str2map=function(str){
    input=unlist(strsplit(str, '[[:space:]]+'))
    chrname=input[1]
    input=input[-(1:2)]
    markers=input[seq(1,(length(input)-1), by=2)]
    positions=as.numeric(input[seq(2,length(input), by=2)])
    #output=cbind(markers, positions)
    output=positions
    names(output)=markers
    #colnames(output)=c("Marker", paste("Chr", chrname, "(cM)"))
    attr(output, "Chr")=chrname
    class(output)=c("A", "map")
    return(output)
  }

  assignattrs=function(x, attrs) {
    for (param in names(attrs)) {
      if(length(attrs[[param]])<1) next
     # print(paste("Param name:", param))
     # print(paste("Param value:", attrs[[param]]))
      attr(x, param)=attrs[[param]]
    }
    return(x)
  } # assignattrs
  
qtlcart2mscanone=function(filenames, traits, methods, offset=0, column="H0:H3", lod.transform=TRUE) {
  output=list()
  ntraits=length(filenames)
  for (i in 1:ntraits) {
    qtlobj=import.qtlcart(filenames[i], methods[i])
    scanone_obj=qtlobj[,c(1,3, match(column, colnames(qtlobj)))]
    colnames(scanone_obj)=c("chr", "pos", "lod")
    scanone_obj[,2]=scanone_obj[,2]*100 + offset # convert to cM
    if(lod.transform) scanone_obj[,3]=scanone_obj[,3]/4.61
    attr(scanone_obj, "model")="normal"
    attr(scanone_obj, "type")="f2"
    attr(scanone_obj, "method")=attr(qtlobj, "method")
    output[[i]]=scanone_obj
    # Asssuming that all traits have the same map
    if(i==ntraits) {
      map=attr(qtlobj, "map")
      if(length(offset)==1) offset=rep(offset, length(map))
      for(i in length(map)) {
        map[[i]]=map[[i]]+offset[i]
      }
      attr(output, "map")=map
    } # if last trait
  } # for i
  attr(output, "method")=methods[1]
  class(output)="mscanone"
  names(output)=traits
  return(output)
} # qtlcart2mscanone

# Plot Scanone restuls in MScanone Objects  V 1.1
# Description
# Plot fucntion for mscanone object, which are created with the trai.scan funtion and hold multiple scanone results
# Usage
# plot.qtl(mscanone, plotname="Qtlplot", method=attr(mscanone, "method"), phenotypes=names(mscanone)
#               title=plotname, addname=attr(mscanone, "addname"), intname=attr(mscanone, "intname"),
#               lodint=TRUE, lodintfun="bayesint", lodtype="lines", mtitle=TRUE, mtick="line", incl.markers=TRUE, incl.markers.names=TRUE,
#               which.markers=NULL, lwd=2, xlab="Map position (cM)", overlay=T, device=getOption("device"), width=700,
#               height=700, xlim=NULL, ylim=NULL, las=1, markers.cex=.9, cex=par("cex"), device.scale=1, legvpos=0.26, cutoff=0, ...)
#                
# Arguments
#       mscanone     			Object of class mscanone, a list of scanone results
#       method       			method used in scanone. Optional. It is taken from the 'method' attribute of mscanone object
#		phenotype				string or vector of phenotype names to plot
#       plotname     			character string used to name the output graphic file
#       title        			character string used to title the fiugure
#       addname      			character vector of additive covariate names
#       intname      			character vector of interacting covariate names
#       mtitle       			if TRUE, the full name of the method used is added to the title
#       lodint					if TRUE, lod intervals will be plotted
#       lodintfun				name of function to calculate lod intervals
#       lodtype                 type of lodinterval symbol c("lines", "bars")
#       mtick                   type of tickmarks for markers c("line", "triangle")
#       incl.markers            if TRUE, marker tickmarks will be included
#       incl.markers.names		if TRUE, marker names will be included
#       which.markers			vector of indices of markers to include
#       overlay      			if TRUE, eahc trait is plotted on a single graph with different colors
#       lwd,xlab,xlim,ylim, las	customizable graph arguments passed to plot function
#       width, height 			Dimensions of the png graphic device
#       cutoff       			minimum peak LOD to declare QTL and draw a CI line. Can be scalar or vector
#       device.scale			ratio between the size of the active graphic device and the default
#		traits.labels			labels to use for traits
#       ...          furher arguments passed to the plot function
# V 1.1 modified to match legend colors with plot colors  RJG.

plot.mscanone=function(mscanone, plotname="Qtlplot", method=attr(mscanone, "method"), phenotypes=names(mscanone),
               title=plotname, addname=attr(mscanone, "addname"), intname=attr(mscanone, "intname"), incl.legend=TRUE,
               lodint=TRUE, lodintfun="mlodint", lodint.args=list(sep=1), lodtype="lines", col=1:length(phenotypes), mtitle=TRUE, mtick="line", incl.markers=TRUE, incl.markers.names=TRUE,
               which.markers=NULL, lwd=2, xlab="Map position (cM)", overlay=T, device=getOption("device"), outer=FALSE, width=700,
               height=700, xlim=NULL, ylim=NULL, las=1, markers.cex=.9, cex=par("cex"), device.scale=1, legvpos=0.26, leghpos=NULL, cutoff=0, ...)  {
  require("qtl") || stop("qtl package not available.")
  if(length(names(mscanone))==0) stop("First argument in not an mscanone object")
  if(all(is.na(pmatch(phenotypes, names(mscanone))))) stop(paste("Phenotypes not included in mscanone object: ", paste(phenotypes, collapse=", "), ".", sep=""))
  op<-list()
  nphenos=length(phenotypes)
  nphenos<=6 || stop("Maximun number of graps per plot is 6.")
  ncols=1
  nrows=1
  toplns=4.1
  bottlns=10
  traitcolors=col
  if(nphenos>2){
    ncols=2
        nrows=ceiling(nphenos/2)
  }
   
  if(nphenos>1){ncols=2}

  methodname=switch(method,
                        em="Expectation Maximization",
                        mr="Marker Regression",
                        hk="Haley-Knott",
                        "mr-imp"="Single imputation",
                        "mr-argmax"="Viterbi algorithm",
                        "CIM"="Composite Interval Mapping")
 if (mtitle) {
   title=paste(title, methodname, "Scan")
 }
 if(overlay & incl.legend) {
  #get(device)(filename=paste(plotname,".png", sep=""), width=width,  height=height)
   op=c(op,par(mar=c(bottlns,4.1,toplns,2.1)))
 } else { # if not overlay
    #png(paste(plotname,".png", sep=""), width=330*ncols,  height=265*nrows+100)
    if(mtitle){
      if(is.null(addcov)){
        op=c(op,par(mfrow=c(nrows,ncols), oma=c(0,0,2,0)))
      } else {
        op=c(op,par(mfrow=c(nrows,ncols), oma=c(0,0,3,0)))
      }
    }
 }
 DSize=0.045*cex*markers.cex/device.scale# /(par("mar")[1]/op$mar[1])
 for(chr  in unique(as.character(mscanone[[1]][,1]))) { # loop through chromosomes
 # this assumes that all phenotypes were scanned for the same chromosomes. These are
 # taken arbitrarily only from the first pehnotype.
    xmax=xlim[2]
    xmin=xlim[1]
    ymax=ylim[2]
    ymin=ylim[1]

  if(overlay) {
    if(is.null(xmax)) {
      xmax=0
      xmin=100
    }
    if(is.null(ymax)) {
      ymax=0
      ymin=0
    }
    realymin=ymin
    if(!xmax | !ymax) {
      for (i in phenotypes){
        if(is.null(xlim[2])){xmax=max(xmax, mscanone[[i]][,2])}
        if(is.null(xlim[1])){xmin=min(xmin, mscanone[[i]][,2])}
        if(is.null(ylim[2])){ymax=max(ymax, mscanone[[i]][,3])}
      } # for every phenotype
    } # if not xmax or ymax
    if (lodint) {
     lodints=numeric()
     plottrait=logical()
     if(length(cutoff)==1) cutoff=rep(cutoff, nphenos)
     for(chr  in unique(as.character(mscanone[[phenotypes[1]]][,1]))) { # loop through chromosomes
       # this assumes that all phenotypes were scanned for the same chromosomes. These are
       # taken arbitraily only from the first pehnotype.
       for (i in 1:nphenos){ # loop through phenotypes
         pheno=phenotypes[i]
         scanres=mscanone[[phenotypes[i]]][mscanone[[phenotypes[i]]][,1]==chr,]
         traitint=do.call(lodintfun, c(list(scanres,chr), lodint.args))
         #print(class(traitint))
         #print(paste("Pheno:", pheno))
         #print(paste("Cutoff:", cutoff[i]))
         if("mlodint" %in% class(traitint)) {
            aresig=sapply(traitint, function(x) x[2,"lod"]>=cutoff[i])
            #print(aresig)
            plottrait=append(plottrait, aresig)
            #print(traitint)
            traitint=t(sapply(traitint, function(x) x[c(1,nrow(x)),"pos"]))
            if(!is.null(traitint)){
              #print(paste(phenotypes[i],": ", sep=""))
              #print(traitint)
              lodints=rbind(lodints, traitint)
            }
            colidx=match(traitcolors[i], col)
            if(i==1) {
              col=c(rep(col[1], nrow(traitint)), col[2:length(col)])
            } else {
              if(i!=length(phenotypes))
                col=c(col[1:(colidx-1)], rep(traitcolors[i], length(aresig)), col[(colidx+1):length(col)])
              else
                col=c(col[1:(colidx-1)], rep(traitcolors[i], length(aresig)))
            }  
         } else {

           plottrait=append(plottrait, max(traitint[,3])>=cutoff[i])
           traitint=traitint[c(1, nrow(traitint)),2]
           #if(!is.list(traitint)) traitint=list(traitint)
           #traitint=lapply(traitint, function(x) x[c(1,nrow(x)),2])
           lodints=rbind(lodints, traitint)
         }
         
       } # loop through phenotypes
     } # loop through chrmosomes
     #print(lodints)
     #print(plottrait)
     lodints=lodints[plottrait, ]
     if(sum(plottrait)==1) lodints=t(as.matrix(lodints))
     col=col[plottrait]
     plottrait=plottrait[plottrait]
     #print(cbind(lodints, col))
     #print(plottrait)
     #print(lodints)
     #print(col)
     if(sum(plottrait)>0) {
       intclusters=find.intclusters(lodints)
       #print(intclusters)
       if(lodtype=="lines" & is.null(ylim)) { # add some space on top of the graph for od int lines
         ymax=ymax+ymax*.02*(max(sapply(intclusters, length))+1) # THE NUMERIC VALUE CAN BE ADJUSTED
       }
     }
   } # if lodint
  } # if overlay
  
  i=1

  for (pheno in phenotypes){ # loop through phenotypes
    #pheno=phenotypes[i]
      scanres=mscanone[[pheno]][mscanone[[pheno]][,1]==chr,]
      markers=attr(mscanone,"map")[[chr]]
      if(overlay) {
        if(i==1) {
          nletters=max(nchar(names(markers)))

          if(incl.markers & is.null(ylim)) ymin = 0 # ymin-DSize*nletters*(ymax-ymin)

          plot(as.matrix(scanres)[,2:3], type="l", lwd=lwd, xlab=xlab, xlim=c(xmin,xmax), ylim=c(ymin,ymax), col=traitcolors[i], las=las, ...)
          if(incl.markers) {
            pos=as.numeric(markers)
            if (mtick == "line") {
              rug(pos, 0.02, quiet = TRUE)
            } else {
              a <- par("usr")
              points(pos, rep(a[3] + diff(a[3:4]) * 0.04, length(pos)),
              pch = 17, cex = 1.5)
            }
#            if(incl.markers.names) {
 #             if(is.null(which.markers)) which.markers=1:length(markers)
  #            text(pos[which.markers], ymin, labels=names(markers)[which.markers], srt=90, cex=markers.cex, adj=0)
   #         } # if incl.markers.names
          } # if incl.markers
        } else {
          lines(as.matrix(scanres)[,2:3], type="l", lwd=lwd, col=traitcolors[i], ...)
        } # if first pheno
        #if (lodint) {
        #    traitint=lodint(scanres,chr)
        #    traitint=traitint[c(1, nrow(traitint)),2]
        #    lodints=rbind(lodints, traitint)
        #} # if lodint

      } else {
        plottitl=paste(pheno)
        plot(as.matrix(scanres)[,2:3], main=plottitl, type="l", lwd=lwd, xlab=xlab, ...)
    } # if overlay
    i=i+1
  } # for i
 } # for chr
 
  if(lodint & length(plottrait)>0) {
    ### print(lodints) ####
    # plot.int(lodints, col=1:nphenos, ymax=ymax, type=lodtype)
      plot.int(lodints, intclusters, col=col, ymax=ymax, type=lodtype, plottrait)
  }
  if(mtitle){
    if(!is.null(addname)){
      title=paste(title, "\nwith", addname, "addcovar")
    }
    if(!is.null(intname)){
      title=paste(title,"and",intname, "intcovar")
    }
  }
  if(overlay ) {
    title(title)
    if(incl.legend){
      op=c(op,par(xpd=TRUE))
      nlines=12
      ypos=ymin-DSize*nlines*(ymax-realymin)
      if(is.null(leghpos)){
        leghpos=(xmax+xmin)*.5
      }
      if(outer){
        tmp <- cnvrt.coords(.9,.7, 'tdev')$usr
        legend(tmp, phenotypes, lwd=2, col=col, ncol=ceiling(nphenos/2), xjust=0.5, yjust=1)
      } else {
        legend(leghpos, ymin-(ymax-ymin)*legvpos, phenotypes, lwd=2, col=col, ncol=ceiling(nphenos/2), xjust=0.5, yjust=1)
      }
    }
  } else {
    mtext(title, outer=TRUE)
  }
  par(op)
  #dev.off()
} # plot.mscanone

plot.lodint=function(lodint, lty=2, col=1, type=c("lines", "bars"), y=NULL) {
  lowint=numeric()
  highlim=numeric()
  if(length(dim(lodint))>1) {
    lowlim=lodint[1,2]
    highlim=lodint[nrow(lodint),2]
  } else {
    lowlim=lodint[1]
    highlim=lodint[2]
  }

  if(type[1]=="lines") {
    lines(c(lowlim, highlim), rep(y, 2), col=col, lwd=3)
    
  } else {
    if(type[1]=="bars") {
      abline(v=lowlim, lty=lty, col=col)
      abline(v=highlim, lty=lty, col=col)
    }
  }  
}

old.plot.int=function(intervals, col, ymax, type) {
  
  intclusters=find.intclusters(intervals)
  for(i in 1:length(intclusters)){
    plot.lodclust(intervals[intclusters[[i]],], col=col[intclusters[[i]]], type, ymax)
  }
} # old.plot.int

plot.int=function(intervals, intclusters, col, ymax, type, plotlog) {
  for(i in 1:length(intclusters)){
    incluster=intclusters[[i]]
    logincluster=plotlog[incluster]
    outidx=incluster[logincluster]
    if(length(outidx)>0)
      plot.lodclust(intervals[outidx,], col=col[outidx], type, ymax)
  }
} # plot.int


plot.lodclust=function(intmat, col, type, ymax) {
  intmat=matrix(intmat, ncol=2)
  onestep=ymax*.02
  y=ymax
#  print(intmat)
  for(i in 1:nrow(intmat)) {
#    print(y)
    plot.lodint(intmat[i,], col=col[i], type=type, y=y)
    y=y-onestep
  }
} # plot.lodclust

find.intclusters=function(intmat){
  overlaps=function(int1, int2) {
#  print(paste("int1[1] - int2", sum(abs(int1[1] - int2)), "int2[2]-int2[1]", (int2[2]-int2[1]) )) ###
#  print(paste("int2[1] - int1", sum(abs(int2[1] - int1)), "int1[2]-int1[1]", (int1[2]-int1[1]))) ###
#    if(sum(abs(int1[1] - int2)) == max((int1[2]-int1[1]), (int2[2]-int2[1])) | sum(abs(int1[2] - int2)) == max((int2[2]-int2[1]), (int2[2]-int2[1])) ) {
    if(abs(sum(abs(int1[1] - int2)) - (int2[2]-int2[1]))<1e-3 | abs(sum(abs(int2[1] - int1)) - (int1[2]-int1[1]))<1e-3 
       | any(abs(int1[2:1]-int2)<1)) {
#      print("  Overlap!!")
      return(TRUE)
    }
  return(FALSE)
  }
  indices=1:nrow(intmat)
  clusters=list()
  i=0
  while(length(indices)>0) {
    i=i+1
#    print(paste("Int1 line:", i))
    members=indices[1]
    secindx=indices[-1]
    for (j in secindx) {
#      print(paste("  int2 line:", j))
      if(overlaps(intmat[indices[1],], intmat[j,])) {
        members=unique(append(members, j))
        indices=indices[indices!=j]
      }
    }
    indices=indices[-1]
    members=members[order(intmat[members,2]-intmat[members,1], decreasing=TRUE)]
    clusters[[i]]=members
  }
  return(clusters)
} # find.intclusters

read_qtl_data=function(input="rqtl_input.csv") {
  require("qtl") || stop("qtl package not available.")
  cross=read.cross(format="csv", file=input, na.strings=c("NA", " ", "-", "?"), genotypes=c("B","H","C","Not B", "Not C"))
  cross=est.rf(cross)
  cross=drop.nullmarkers(cross)
  cross=est.rf(cross)
  orders=ripple(cross, "17", method="countxo")
  bestorder=orders[orders[,ncol(orders)]==min(orders[,ncol(orders)]),]
  
  newmap=est.map(cross, error.prob=.01, verbose=F)
  cross=replace.map(cross, newmap)
  cross=calc.genoprob(cross, error.prob=.01)
  cross=calc.errorlod(cross, error.prob=.01)
  cross=calc.genoprob(cross, step=2, error.prob=.01)
  return(cross)
} # read_qtl_data

qtlplots=function(cross, plotname="Qtlplot", phenotypes, title=plotname, addcov=NULL, addname=NULL, 
                  intname=NULL, intcov=NULL, method="hk", mtitle=TRUE, lwd=2, 
                  xlab="Map position (cM)", add=F, ...)  {
  require("qtl") || stop("qtl package not available.")
  nphenos=length(phenotypes)
  nphenos<=6 || stop("Maximun number of graps per plot is 6.")
  ncols=1
  nrows=1
  if(nphenos>2){
    ncols=2
	nrows=ceiling(nphenos/2)
  }
   
  if(nphenos>1){ncols=2}

  methodname=switch(method,
			em="Expectation Maximization",
			mr="Marker Regression",
			hk="Haley-Knott",
			"mr-imp"="Single imputation",
			"mr-argmax"="Viterbi algorithm")

 if (mtitle) {
   title=paste(title, methodname, "Scan")
 }
  png(paste(plotname,".png", sep=""), width=330*ncols,  height=265*nrows+100)
  if(is.null(addcov)){
    par(mfrow=c(nrows,ncols), oma=c(0,0,2,0))
  } else {
    par(mfrow=c(nrows,ncols), oma=c(0,0,3,0))
  }
  
  for (i in 1:nphenos){
    pheno=phenotypes[i]
    scanres=scanone(cross, pheno.col=find.pheno(cross,pheno), method=method,
             addcovar=addcov, intcovar=intcov )
    #plottitl=paste(pheno, methodname, "Scan")
    plottitl=paste(pheno)
    plot(as.matrix(scanres)[,2:3], main=plottitl, type="l", lwd=lwd, xlab=xlab, add=add, ...)
    varadd=add-1
  }
  if(!is.null(addcov)){
    title=paste(title, "\nwith", addname, "addcovar")
   }
    if(!is.null(intcov)){
    title=paste(title,"and",intname, "intcovar")
   }

  mtext(title, outer=TRUE)                                    
  dev.off()
}

# qtlint applies a test to discrimiate the most likely side from a middle
# marker where the QTL is located
qtlint=function(cross, mid, pheno, chr) {
  "cross" %in% class(cross) || stop("Argument cross should be of class cross (see scanone).")
  require(qtl) || stop("R/QTL package not vailable.")
  genotypes=pull.geno(cross, chr)
  #midpos=find.marker.pos(cross, chr, mid)
  flanking=find.flanking.markers(cross, chr, mid)
  
} # qtlint

trait.scan=function(cross, phenotypes, addcov=NULL, addname=NULL, 
                  intname=NULL, intcov=NULL, method="hk", ...) {
                  
  require("qtl") || stop("qtl package not available.")
  if(is.null(intname) & !is.null(colnames(intcov))) intname=colnames(intcov)
  if(is.null(addname) & !is.null(colnames(addcov))) addname=colnames(addcov)

  nphenos=length(phenotypes)
  output=list(0)
         
  for (i in 1:nphenos){
    pheno=phenotypes[i]
    output[[i]]=scanone(cross, pheno.col=find.pheno(cross,pheno), method=method,
             addcovar=addcov, intcovar=intcov, ... )
    names(output)[length(output)]=pheno
    }
  class(output)="mscanone"
  attr(output, "method")=method
  attr(output, "addname")=addname
  attr(output, "intname")=intname
  attr(output, "map")=pull.map(cross)
  
  return(output)
  } # trait.scan

mcim=function(cross, phenotypes, method="hk", cim.fun=cim, ...) {
                  
  require("qtl") || stop("qtl package not available.")
  nphenos=length(phenotypes)
  output=list(0)
         
  for (i in 1:nphenos){
    pheno=phenotypes[i]
    output[[i]]=cim.fun(cross, pheno.col=find.pheno(cross,pheno), method=method, ... )
    names(output)[length(output)]=pheno
    }
  class(output)="mscanone"
  attr(output, "method")=method
  attr(output, "map")=pull.map(cross)

  return(output)
  } # mcim


bayesint.mscanone=function(x, chr, prob=0.95, cutoff=0, collapse=T, decimals=NULL) {
  # x	object of  class mscanone
  # chr	chromosome
  require(qtl) || stop("qtl package not available.")
  mint=lapply(x, bayesint, chr, prob)
  if(length(cutoff)==1) cutoff=rep(cutoff, length(x))
  out=numeric()
  for(i in 1:length(mint)) {
    if(max(mint[[i]][,3])>=cutoff[i]){
      rowout=as.vector(lodint2lims(mint[[i]], collapse=collapse, decimals=decimals))
      out=rbind(out, rowout)
      rownames(out)[nrow(out)]=names(x)[i]
    }
  }
    #out=t(sapply(mint[[, lodint2lims))
  if(collapse) colnames(out)=paste(prob*100, "%", " CI", sep="")
  return(out)
} # bayesint.mscanone

lodint2lims=function(lodint, cutoff=0, collapse=F, decimals=NULL, with.peak=F, with.lod=F){
  if(length(dim(lodint))>1) {
    if(max(lodint[,3]>=cutoff)){
      lowlim=lodint[1,2]
      highlim=lodint[nrow(lodint),2]
      peak=lodint[2, "pos"]
      lod=lodint[2, "lod"]
      chr=round(as.numeric(as.character(lodint[2, "chr"])), 0)
    } else {
      return(NULL)
    }
  } else {
    lowlim=lodint[1]
    highlim=lodint[2]
  }
  out=c(lowlim, highlim)
  
  names(out)=c("CI.low", "CI.up")
  if(!is.null(decimals)) out=format(round(out, decimals), digits=decimals+3)
  if(with.lod) {
    out=append(lod, out)
    names(out)[1]="lod"
  }
  if(with.peak) {
    out=append(c(chr, peak), out)
    names(out)[1:2]=c("chr", "peak")
  }
  class(out)=c("numeric", "lims")
  if(collapse) out=paste(out, collapse=" - ")
  return(out)
} # lodint2lims

mlodint= function(x, ...) 
			UseMethod("mlodint")

mlodint.scanone=function(scanone, chr, drop=1.5, sep=5, cutoff=0, conservative=T) {
# Finds multiple lod intervals in one chromosome
#	scanone		object of clas scanone
#	chr			chromosome name
#	drop		difference in LOD between peak and valley
#	sep			minimum separation between peaks
#	cutoff		LOD threshold to declare a QTL
  peaks=mpeaks(scanone, sep, cutoff)
  peaks=peaks[peaks[,"chr"]==chr,]
 #print(peaks)
  if(length(peaks)==0) return(NULL)
 #print(peaks)
  valley=peaks[,"lod"]-drop
 #print(valley)
  #sigregion=sapply(scanone[,"lod"],'>=',valley)
  sigregion=matrix(rep(scanone[,"lod"], length(valley)), nrow=length(valley), byrow=T)
 # print(sigregion)
  if(conservative)
    sigregion=(sigregion-valley)>0
  else
    sigregion=(sigregion-valley)>=0
 #print(dim(scanone))
 #print(colnames(scanone))
 #print(scanone[,"lod"])
 #print(sigregion)
  change=abs(t(diff(t(sigregion))))
 #print(dim(change))
  change=cbind(rep(0, nrow(peaks)), change)
 #print(change)
 maxidx=match(peaks[,"pos"], scanone[,"pos"])
 #print(maxidx)
 updiffs=list()
 downdiffs=list()
 contregions=matrix(nrow=dim(change)[1], ncol=dim(change)[2])
 #marklims=c(1,nrow(scanone))
 for(i in 1:length(maxidx)){
   downdiffs=logical()
   updiffs=cumsum(change[i,(maxidx[i]+1):nrow(scanone)])==0
   if(maxidx[i]!=1) { # peak not in begining of chr
     #downdiffs=cumsum(change[i, maxidx[i]:1])==0
     downdiffs=cumsum(change[i, maxidx[i]:2])==0
     downdiffs=downdiffs[length(downdiffs):1]
     #downdiffs=downdiffs[-length(downdiffs)]
     #controw=c(downdiffs[-length(downdiffs)], rep(TRUE,3), updiffs[-1]) # most inclusive interval
     if(conservative)
       controw=c(downdiffs[-1], rep(TRUE,3), updiffs[-length(updiffs)]) # most inclusive interval
     else       
       controw=c(downdiffs, TRUE, updiffs)
     contregions[i,]=controw
   } else {
     controw=c(TRUE, updiffs)
     contregions[i,]=controw
   }
 }
 #print(contregions)
 intregions=apply(contregions, 1, function(x) scanone[x,"pos"])
 if(length(maxidx)==1) intregions=list(as.vector(intregions)) # security check for special beaviour of apply
 #print(intregions)
 if("matrix" %in% class(intregions))
   intervals=t(apply(intregions, 2, range))
 else
   intervals=t(sapply(intregions, range))
 #print(intervals)
 lodints=ints2lodints(scanone, peaks, intervals, chr)
 class(lodints)=c("mlodint", "list")
 #print(lodints)
 lodints=unique(lodints)
 return(lodints)
} # mlodint

unique.mlodint=function(mlodint){
  inidx=1:length(mlodint)
  repidx=numeric()
  if(length(inidx)>1){ 
    for (i in inidx){
      other=inidx[-c(i, repidx)]
      for(j in other){
        int1=lodint2lims(mlodint[[i]])
        int2=lodint2lims(mlodint[[j]])
        # Lets find if this interval is contained within any other one or if they are exacly equal
        if(all(int1==int2) | (abs(sum(abs(int1[1] - int2)) - (int2[2]-int2[1]))<1e-3 & abs(sum(abs(int1[2] - int2)) - (int2[2]-int2[1]))<1e-3 )) {
        #if(all(mlodint[[i]][c(1,3),"pos"]==mlodint[[j]][c(1,3),"pos"])) {
          #print(paste("rep:", i))
          # if true, save the one with highest lod score
          if(mlodint[[i]][2,"lod"]>mlodint[[j]][c(2),"lod"])
            repidx=append(repidx, j)
          else
            repidx=append(repidx, i)
        }
      }
    }
  } else {
    return(mlodint)
  }
  if(length(repidx)==0) return(mlodint)
  output=mlodint[-repidx]
  if(length(output)==0) return(NULL)
  class(output)=c("mlodint", "list")
  return(output)
} # unique.mlodint

mlodint.mscanone=function(mscanone, chr, drop=1.5, sep=5, cutoff=0, collapse=F, decimals=NULL, with.peak=T, with.lod=T) {
#print("Fx: mlodint.mscaone")
  # x	object of  class mscanone
  # chr	chromosome
  if(length(cutoff)==1) cutoff=rep(cutoff, length(mscanone))
  mint=list()
  for(i in 1:length(mscanone)){
    mint[[i]]=mlodint(mscanone[[i]], chr, drop, sep, cutoff[i])
  }
  #mint=lapply(mscanone, mlodint, chr, drop, sep, cutoff)
  out=numeric()
  rowlabs=character()
  if(length(mint)>0) { # if any interval
    for(i in 1:length(mint)) {# for each trait
      if(length(mint[[i]])>0){ # if any peaks
        for(j in 1:length(mint[[i]])) { # for each peak
        #print(paste(names(mscanone)[i], ": peak", j))
          if(max(mint[[i]][[j]][,3])>=cutoff[i]){
            rowout=lodint2lims(mint[[i]][[j]], collapse=collapse, decimals=decimals, with.peak=with.peak, with.lod=with.lod)
            out=rbind(out, rowout)
            rowlabs=c(rowlabs, names(mscanone)[i])
          }  
        } 
      }
    }
  } else {
    print("No significant QTL found.")
    return(NULL)
  } 

  #out=t(sapply(mint[[, lodint2lims))
  rownames(out)=dimlabels(rowlabs)
#  colnames(out)=names(rowout)
  if(collapse) colnames(out)="CI"
  return(out)
} # mlodint.mscanone

ints2lodints=function(scanone, qtl, ints, chr) {
# Interal function that conversts a matrix of intervals dim(qtl,2) into a list of lodint objects
   out=list()
   for(i in 1:nrow(qtl))
     out[[i]]=int2lodint(ints[i,], qtl[i, "pos"], scanone, chr)
   return(out)
} # ints2lodints

int2lodint=function(int, peak, scanone, chr) {
   pos=c(int[1], peak, int[2])
   lod=scanone[scanone[,"chr"]==chr,]
   lod=lod[nummatch(pos,lod[,"pos"], .0001),"lod"]
   out=data.frame(chr=rep(chr, 3), pos, lod)
   markers=rownames(scanone)[nummatch(pos, scanone[,"pos"], .0001)]
   rownames(out)=dimlabels(markers)
   return(out)
} # int2lodint

dimlabels=function(x) {
   out=x
   for (i in 2:length(x)) {
     prev=x[1:(i-1)]
     if(x[i]%in%prev) out[i]=paste(x[i], sum(x[i]==prev), sep=".")
   }
  return(out)
} #dimlabels
nummatch=function(x, table, precision=1e-10){
   out=numeric()
   for(i in 1:length(x)) out[i]=(1:length(table))[abs(x[i]-table)<precision][1]
   return(out)
} # nummatch

mpeaks= function(x, ...) 
			UseMethod("mpeaks")

mpeaks.scanone=function(scanone, sep=5, cutoff=0) {
# Internal function to find multple peaks
   chrout=character()
   output.lod=numeric()
   output.pos=numeric()
   qtllabs=character()
 for (c in unique(as.character(scanone[,"chr"]))) { # for each chr
   #print(paste("chr:", c))
   out=numeric()
   outlod=numeric()
   max=-1
   prev=0
   prevmax=0
   chrscan=scanone[scanone[,"chr"]==c,]
   for(i in 1:nrow(chrscan)) {
     markers=rownames(scanone)
     pos=chrscan[i,"pos"]
     lod=chrscan[i,"lod"]
     #print(paste("pos:", pos, "lod:", lod))
     if(max>=lod){ # if going down or plateau
       if((i < 3) || (prevmax<max) || (max==lod & i!=nrow(chrscan) & chrscan[i+1,"lod"]<=lod)) { # if after marker 2 OR peak OR plateau and valley later
         if(max>=cutoff) {# is significant
           if(length(out)==0 || (prev-out[length(out)])>=sep) { # at good distance
             out=c(out, prev)
             outlod=c(outlod, max)
             qtllabs=c(qtllabs, markers[i-1])
           } else { # too close
             if(max>outlod[length(out)]) { # this peak higher than previous one
               out[length(out)]=prev
               outlod[length(outlod)]=max
               qtllabs[length(qtllabs)]=markers[i-1]
             }
           }
         }
       }
     }
     # Advance one step
     prevmax=max
     max=lod
     prev=pos
   }
   chrout=c(chrout, rep(c, length(outlod)))
   output.lod=c(output.lod, outlod)
   output.pos=c(output.pos, out)
 } # each chr
   #print(out)
   if(length(output.lod)>0) {
     qtllabs=dimlabels(qtllabs)
     output=data.frame(chr=chrout, pos=output.pos, lod=output.lod)
     rownames(output)=qtllabs
     return(output)
   } else NULL
} # mpeaks

mpeaks.mscanone=function(mscanone, chr, sep, cutoff) {
  if(length(cutoff)==1) cutoff=rep(cutoff, length(mscanone))
  #out=matrix(nrow=1,ncol=3)
  #colnames(out)=c("chr", "pos", "lod")
  out=character()
  outrowlabs=character()
  for (i in 1:length(mscanone)){
    tmp=mpeaks(mscanone[[i]], sep, cutoff[i])
    if(is.null(tmp)) next
    outrowlabs=c(outrowlabs, rep(names(mscanone)[i], nrow(tmp)))
    out=rbind(out, tmp)
  }
  outrowlabs=dimlabels(outrowlabs)
  rownames(out)=outrowlabs
  return(out)
}

summary.mscanone=function(x, cutoff) {
# x			object of  class mscanone
# cutoff	minimum peak LOD to declare QTL and draw a CI line. Can be scalar or vector
require(qtl) || stop("qtl package not available.")
out=numeric()
if(length(cutoff)==1) cutoff=rep(cutoff, length(x))
for (i in 1:length(x)){
  qtllist=summary(x[[i]], cutoff[i])
  #print(qtllist)
  if(any(!is.na(qtllist))){
    out=rbind(out, qtllist)
    rownames(out)[nrow(out)]=names(x)[i]
  }
}
if(length(dim(out)>0)){
  colnames(out)=c("chr", "pos", "lod")
}
return(out)
} # summary.mscanone

lims2flank=function(x, map, chr, collapse=TRUE) {
#  x	a lims object. Tipically the output of lodint2lims
   if(dim(x)[1]==1){
     x=t(as.matrix(x[,c("CI.low", "CI.up")]))
     indx=find.flanks(x, map[[chr]], collapse)
   } else 
       indx=apply(x[,c("CI.low", "CI.up")], 1, find.flanks, map[[chr]], collapse)
   if(!collapse) indx=t(indx)
   return(indx)
} # lims2outflank

find.flanks=function(lims, mapvector, collapse=F) {
  if(length(lims)==1) {
    if(grep("-", lims)) lims=as.numeric(unlist(strsplit(lims, "[[:space:]]*-[[:space:]]*")))
    if(length(lims)!=2) stop("argument lims is not a limits object")
  }
  mapvector=sort(mapvector)
  low=(1:length(mapvector))[mapvector<=lims[1]]
  low=low[length(low)]
  up=(1:length(mapvector))[mapvector>=lims[2]]
  up=up[1]
  #return(c(low, up))
  out=names(mapvector)[c(low, up)]
  names(out)=c("flank.low", "flank.up")
  if(collapse) out=paste(out[1], out[2], sep=" - ") 
  return(out)
} # find.flanks

mfitqtl=function(cross, phenos, qtl, covar=NULL, formula, method=c("imp"), dropone=TRUE, 
        get.ests=FALSE, stats=c("df", "SS", "MS", "LOD", "%var", "Pvalue(Chi2)", "Pvalue(F)")) {
#print("Fx: mfitqtl")
# phenos	a vector of phenotype names
# qtl 		An object of class qtl (or a vector of), as output from makeqtl.
# covar 	A data.frame of covariates
# formula	An object of class formula (or a vector of) indicating the model to be fitted. 
#           QTLs are referred to as Q1, Q2, etc. Covariates are referred to by their names
#           in the data frame covar.
# method 	Indicates whether to use the EM algorithm or imputation. (Only imputation is
#			implemented at this point.)
# dropone 	If TRUE, do drop-one-term analysis.
# get.ests 	If TRUE, return estimated QTL effects and their estimated variance-covariance
#			matrix.
# stats		a vector of statistics names to include in the report

  require(qtl) || stop("qtl package not available.")
  if(length(qtl)==1) qtl=as.list(rep(qtl, length(phenos)))
  indeps=attr(terms(formula), "term.labels")
  out=numeric()
  outcols=character()
  for (i in 1:length(phenos)) {
    #print(paste(i,":", phenos[i]))
    y=cross$pheno[,find.pheno(cross, phenos[i])]
    fitqtl.obj=fitqtl(y, qtl[[i]], covar, formula, method, dropone, get.ests)
    #print(fitqtl.obj)
    if(dropone & "result.drop" %in% names(fitqtl.obj)) {
      #print(str(fitqtl.obj[["result.drop"]]))
      if("SS" %in% stats)
        stats[grep("SS", stats)]="Type III SS"
      if("MS"%in% stats)
        stats=stats[-grep("MS", stats)]
      outrow=matrix(t(fitqtl.obj[["result.drop"]][,stats]), nrow=1)
      #print(rownames(fitqtl.obj[["result.drop"]]))
      outcols=sapply(rownames(fitqtl.obj[["result.drop"]]), paste, stats, sep=".") 
    } else {
      #print(stats)
      #print(colnames(fitqtl.obj[["result.full"]]))
      #print(rownames(fitqtl.obj[["result.full"]]))
      outrow=matrix(t(fitqtl.obj[["result.full"]]["Model",stats]), nrow=1)
      outcols=stats
    }
    if(get.ests) {
   #print("fitqtl will provide eststimates")
      tmp=summary(fitqtl.obj)
      outrow=append(outrow, fitqtl.obj$ests$ests)
      estnames=names(fitqtl.obj$ests$ests)
      ests.se=sqrt(fitqtl.obj$ests$covar[diag(nrow(fitqtl.obj$ests$covar))>0])
      estscols=c(estnames, paste(estnames, "se", sep="."))
      outrow=append(outrow, ests.se)
    }
    out=rbind(out, outrow)
  } # for pheno
  
  #outcols=as.vector(sapply(indeps, paste, stats, sep="."))
  if(get.ests) outcols=append(outcols, estscols)
  
#  out=as.data.frame(out)
  rownames(out)=dimlabels(phenos)
  colnames(out)=outcols
  return(out) 
} # mfitqtl

mmakeqtl=function(cross, chr, pos, qtl.name=NULL) {
# cross		An object of class cross. See read.cross for details.
# chr 		Vector indicating the chromosome for each QTL.
# pos 		Vector (of same length as chr) indicating the positions on the chromosome to
#			be taken. If there�s no marker or pseudomarker at a position, genotypes for the
#			nearest positions are taken.
# qtl.name 	The user-specified name for each QTL, used in the drop-one-term ANOVA table
#			in fitqtl. If unspecified, the names will be of the form "Chr1@10" for a
#			QTL on Chromsome 1 at 10 cM.
  require(qtl) || stop("qtl package not available.")
  if(length(chr)<length(pos)) chr=as.list(rep(chr, length(pos)))
  if(length(pos)<length(chr)) pos=as.list(rep(pos, length(chr)))
  if(length(qtl.name)==1) qtl.name=as.list(rep(qtl.name, length(chr)))
  out=list()
  for(i in 1:length(pos)){
    #rint(chr[[i]])
    #print(pos[[i]])
    if(is.null(qtl.name))
      qtl.obj=makeqtl(cross, chr[[i]], pos[[i]])
    else 
      qtl.obj=makeqtl(cross, chr[[i]], pos[[i]], qtl.name[[i]])
    #print(qtl.obj)
    out[[i]]=qtl.obj
    class(out[[i]])="qtl"
  } # for i
  return(out)
} # mmakeqtl

scanbind=function(..., map=NULL, method=NULL){
  input=list(...)
  args=unlist(input, recursive=F)
  class(args)=c("mscanone", "list")
  if(is.null(map) & !is.null(attr(input[[1]], "map")))
     attr(args, "map")=attr(input[[1]], "map")
  if(is.null(method) & !is.null(attr(args[[1]], "method")))
    attr(args, "method")=attr(args[[1]], "method")
     
  if(!is.null(map))
     attr(args, "map")=map
  if(!is.null(method))
     attr(args, "method")=method
  return(args)
}
is.sig=function(mscanone, cutoff, sep=1) {
  out=logical()
  for (i in 1:length(mscanone)){
    #out=c(out, max(mscanone[[i]][,3])>cutoff[i])
    peaks=mpeaks(mscanone[[i]], sep=sep, cutoff=cutoff[i])
    if(!is.null(peaks))
      out=c(out, rep(i, nrow(peaks)))
  }
  return(out)
}

qtleffects.cross=function(cross, formula, phenotype, position, chromosome, addcovar=NULL, intcovar=NULL, 
                          addnames=NULL, intnames=NULL, stats=c("Estimate","Std. Error","t value","Pr(>|t|)")) {
#print("Fx: qtleffects.cross")                          
#		cross		cross object
#		phenotype	character or vector of phenotype names
#		position	position (in same unit as map in cross) to test. Can be a vector.
#		chromosme   character string or vector with names of chromosome for each position.
#		addcovar	matrix of additive covariables
#		incovars	matrix of intrating covariables
#		addnames	vector of names of additive covariables. Default: colnames(addcovars)
#		intnames	vector of names of interacting covariables. Default: colnames(intcovars)


# Variables are splited by combiantions of levels in int covars and a linear model is fitted. FUN is applied
# to the resutls and a matrix of dim((parameters returned by FUN x interaciton levesl) x phenotypes) is returned.
require(qtl) || stop ("qtl library not available.")
# Security Checks
  if(is.null(addcovar) & !is.null(addnames)){
    for(var in addnames) if(!var%in%colnames(cross$pheno)) stop(paste("addcovar", var,"not found in cross."))
    addcovar=data.frame(cross$pheno[,addnames])
  } else {
    if(is.null(addnames) & !is.null(colnames(addcovar))) addnames=colnames(addcovar)
      #else stop("addcovar or addnames must by provided")
  }
  if(is.null(intcovar) & !is.null(intnames)) {
    for(var in intnames) if(!var%in%colnames(cross$pheno)) stop(paste("intcovar", var,"not found in cross."))
    intcovar=data.frame(cross$pheno[,intnames])
  } else {
    if(is.null(intnames) & !is.null(colnames(intcovar))) intnames=colnames(intcovar)
      #else stop("intcovar or intnames must be provided")
  }  
  for(var in phenotype) if(!var%in%colnames(cross$pheno)) stop(paste("Phenotype",var, "not found."))
  # Recycle variables
  #  to be done...
  # Eliminate additive covars present in intcovars
  # Eliminate additive terms from formula that are present intcovar
  # Eliminate interaction terms from formula that are present in intcovar
  # This assumes that intcovar matrix correctly includes all interacting variables in formula
  elim=colnames(addcovar) %in% colnames(intcovar)
#print(formula)
  formulavec=as.character(formula)
  formulaterms=terms(formula)
  indepterms=attr(formulaterms, "term.labels")
  addinint=indepterms %in% colnames(intcovar)
  inputaddcovars=colnames(addcovar)
#print(paste("intcovar:", paste(colnames(intcovar), collapse=", ")))
#print(paste("addcovar:", paste(colnames(addcovar), collapse=", ")))
#print(indepterms)
  intinint=sapply(colnames(intcovar), regexpr, indepterms)>-1
  if(length(dim(intinint))==0)
    intinint=t(as.matrix(intinint))
#print(intinint)
  intinint=apply(intinint, 1, any)
  if(sum(addinint | intinint)<length(indepterms))
     formula=as.formula(paste("y", paste(indepterms[!addinint & !intinint], collapse=" + "), sep=" ~ "))
  else
     formula=NULL
#print(formula)  
#print(paste("Input addcovars:", paste(colnames(addcovar), collapse=", ")))
  if(sum(elim)<length(elim)) {
    addcovar=as.matrix(addcovar[,!elim])
    if(length(elim)-sum(elim)==1)
       colnames(addcovar)=inputaddcovars[!elim]
#print(paste("addcovars left:", paste(colnames(addcovar), collapse=", ")))
    }
  else
    addcovar=NULL
  nphenos=length(phenotype)
  if(nphenos!=length(chromosome) | nphenos!=length(position)) stop("Arguments chromosome, position, and phenotype must be of same length.")
  if(!is.null(intcovar))
    levellabs=factorints(intcovar, intnames)
  fchar <- paste("pheno ~",paste(c(as.character(formula)[length(formula)], "a + d"), collapse=" + "))
#print(fchar)
# Loop through phenotypes  
 # Find marker close to pos
  out=numeric()
  for(i in 1:nphenos) {
    pheno=cross$pheno[,phenotype[i]]
    marker=pos2marker(position[i], chromosome[i], pull.map(cross))
    a=cross$geno[[chromosome[i]]]$data[,marker]
    d=as.numeric(a==2)
    phendata=data.frame(pheno, a, d)
    if(!is.null(addcovar))
      phendata=cbind(phendata, addcovar)
    if(length(addnames)==1) colnames(phendata)[4]=addnames
    
#print(paste("Final data set vars:", paste(colnames(phendata), collapse=", ")))
#print(paste("Final data set dim:", paste(dim(phendata), collapse=", ")))
  #print(summary(phendata))    
    tmpformula=as.formula(fchar)
    if(!is.null(intcovar)){
      intcovar=as.data.frame(intcovar)    
      tmp <- by(phendata, intcovar, function(x) summary(lm(tmpformula, data=x)))
      }
    else
      tmp=list(summary(lm(tmpformula, data=phendata)))
  #print(tmp)      
    ests=as.vector(sapply(tmp, function(x) x$coeff[c("a","d"),stats]))
    out=rbind(out,ests)
  }
  estnames=paste(rep(stats, each=2), c("[a]","[d]"), sep="") 
  if(!is.null(intcovar))
    estnames=paste(rep(levellabs, each=length(estnames)), rep(estnames, length(levellabs)), sep="_")
  out=as.data.frame(out)
  outrownames=dimlabels(phenotype)
  rownames(out)=outrownames
  colnames(out)=estnames
  return(out)
} # qtleffects.cross

pos2marker= function(x, ...) 
			UseMethod("pos2marker")


pos2marker.default=function(pos, chr, map) {
  diffs=sapply(pos, function(x) abs(x-map[[chr]]))
  if(!is.null(dim(diffs))) # is matrix
    minidx=apply(diffs, 2, function(x) match(min(x), x))
  else
    minidx=match(min(diffs), diffs)
  return(names(map[[chr]])[minidx])
} # pos2marker

pos2marker.matrix=function(pos, chr, map) {
  apply(pos, 2, pos2marker.default, chr, map)
} # pos2marker


old.pos2marker=function(pos, chr, map) {
  diffs=abs(pos-map[[chr]])
  minidx=match(min(diffs), diffs)
  return(names(map[[chr]])[minidx])
} # pos2marker

                      
factorints=function(factors, names=colnames(factors)) {
  factors=as.matrix(factors)
  levellabs=list()
  for(i in 1:length(names)) { 
    validlevs=as.factor(as.character(factors[!is.na(factors[,i]),i]))
    levellabs[[i]]=paste(names[i], levels(validlevs), sep="")
  }
  tmp=character()
  if(length(levellabs)>1) {
    tmp=as.vector(t(sapply(levellabs[[length(levellabs)-1]], paste, levellabs[[length(levellabs)]], sep=":")))
    if(length(levellabs)>2){
      for(i in (length(levellabs)-2):1)
        tmp=as.vector(t(sapply(levellabs[[i]], paste, tmp, sep=":")))
    } # if n factors > 2
  } else {
    return(levellabs[[1]])
  }
  return(tmp)
} # factorints

cim2map=function(scanone, chr=1, pos=2) {
  chromosomes=sort(unique(as.character(scanone[,chr])))
  output=list()
  for (chridx in chromosomes) {
    positions=scanone[scanone[,chr]==chridx,pos]
    markers=rownames(scanone[scanone[,chr]==chridx,])
    names(positions)=markers
    loglik=attr(attr(scanone, "map")[["17"]],"loglik")
    attr(positions,"loglik")=loglik
    output[[chridx]]=positions
  }
  return(output)
} # cim2map

makeQTLtable=function(cross, mscanone, chr, cutoff, formula, qtl.name="Q1", addcovar=NULL, intcovar=NULL, sep=1, map=pull.map(mscaone), lodintfun="mlodint", with.effects=FALSE) {
# Note: you have to include the qtl.name argument (e.g. Q1) in the model
#print("Fx: makeQTLtable")
  get.ests=F
  indepvars=attr(terms(formula), "term.labels")
#print(paste("Indepvars:", paste(indepvars, collapse=", ")))
  interaction=any(regexpr(":", indepvars)>-1)
  realvars=as.logical(regexpr("^[[:alnum:]]+$", indepvars)+1) & indepvars!=qtl.name
  realvars=indepvars[realvars]
#print(paste("Real  vars:", paste(realvars , collapse=", ")))  
  intvars=regexpr(":", indepvars)>-1
  qtl.interaction=regexpr(qtl.name, indepvars)>-1 & intvars
  intvars=gsub(paste(":*", qtl.name, ":*", sep=""), "", indepvars[qtl.interaction])
  usefulvars=indepvars[indepvars!=qtl.name]
  usefulvars=unique(unlist(strsplit(usefulvars, ":", fixed=TRUE)))
  usefulvars=usefulvars[usefulvars!=qtl.name]
#print(paste("Useful vars:", paste(usefulvars, collapse=", ")))
#print(paste("Interaction vars:", paste(intvars, collapse=", ")))
  allvars=unique(c(realvars, intvars))
#print(paste("intvars:", paste(intvars, collapse=", ")))
  if(with.effects & !interaction)
    get.ests=T  
  if(length(realvars)>0 & is.null(addcovar)) {
         colidx=match(usefulvars, colnames(cross$pheno))
         if(any(is.na(colidx))) stop("Additive covariable(s) not found among cross phenotypes.")
         addcovar=as.matrix(cross$pheno[,colidx])
         colnames(addcovar)=usefulvars
   }
  if(length(intvars)>0 & is.null(intcovar)) {
         colidx=match(intvars, colnames(cross$pheno))
         if(any(is.na(colidx))) stop("Interacting covariable(s) not found among cross phenotypes.")
         intcovar=as.matrix(cross$pheno[,colidx])
         colnames(intcovar)=intvars
  }
#print(paste("addcovars:", paste(colnames(addcovar), collapse=", ")))      
#print(paste("intcovars:", paste(colnames(intcovar), collapse=", ")))     
  lodfun=get(lodintfun)
  int=lodfun(mscanone, chr=chr, sep=sep, cutoff=cutoff, collapse=FALSE)
  if(length(int)==0) {
    #print("No significant QTL found.")
    return(NULL)
  }
  #print("there is at least one QTL")
#print("Ok1")  
  flank=lims2flank(int, map=map, chr=chr)
  labs=gsub("(\\.[[:digit:]]+$)", "", rownames(int))
  cutoff.sig=cutoff[pmatch(labs, names(cutoff), duplicates.ok=TRUE)]
  names(cutoff.sig)=rownames(int)
  peak.markers=pos2marker(int[,"peak"], chr, map)
  qtl.model=mmakeqtl(cross, chr=chr, pos=as.list(int[,"peak"]), qtl.name=qtl.name)
#print(paste("Formula for mfitqtl:", formula[3]))
  traits.e.effects=mfitqtl(cross, labs, qtl.model, formula=formula,  covar=addcovar, get.ests=get.ests, dropone=TRUE)
#print("mfitqtl ok")  
  QTL.table=cbind(cutoff=cutoff.sig, int, Flanks=flank, Peak.marker=peak.markers, traits.e.effects)
#print("Ok3")
  if(with.effects & any(qtl.interaction)){
#print("Will use qtleffects.cross") 
    formula=as.formula(paste("y", paste(indepvars[indepvars!=qtl.name & !qtl.interaction], collapse="+"), sep="~"))
    if(length(intvars)>0) {
      intcovar=as.matrix(intcovar[,intvars])
      colnames(intcovar)=intvars                             
    } else
      incovar=NULL
#print(paste("intcovars:", paste(colnames(intcovar), collapse=", ")))
#print(paste("addcovars:", paste(colnames(addcovar), collapse=", ")))
    addvars=colnames(addcovar)
    addcovar=as.matrix(addcovar[,!addvars %in% intvars])
    colnames(addcovar)=addvars[!addvars %in% intvars]
    approxeffects=qtleffects.cross(cross, formula, labs, int[,"peak"], rep(chr, nrow(int)), addcovar=addcovar, intcovar=intcovar)
#print("Ok4")    
    QTL.table=cbind(QTL.table, approxeffects)
  }
#print("Done")  
  return(QTL.table)
}

impcim=function(..., n.imp=50, FUN=mean) {
  imputed=numeric()
  covars=numeric()
  for(i in 1:n.imp) {
    cimres=cim(...)
    if(i==1){ 
      coord=cimres[,1:2]
      covarchr=data.frame(chr=attr(cimres, "marker.covar.pos")[,1])
      covars=attr(cimres, "marker.covar.pos")[,2]
      imputed=cimres[,3]
      next
    }
    imputed=cbind(imputed, cimres[,3])
    covars=cbind(covars, attr(cimres, "marker.covar.pos")[,2])
  }
  #colnames(covars)=NULL
  #colnames(imputed)==NULL

  rownames(covarchr)=paste("covar", 1:nrow(covars), sep="")
  out=apply(imputed, 1, FUN)
  out=cbind(coord, lod=out)
  attr(out, "marker.covar.imp")=covars
  covarpos=cbind(covarchr, apply(covars, 1, FUN))
  attr(out, "marker.covar.pos")=covarpos
  attr(out, "lod.imp")=imputed
  attr(out, "FUN")=as.character(substitute(FUN))
  class(out)=c("scanone", "data.frame")
  return(out)
} # impcim

update.imp=function(imp.scanone,  FUN) {
  output=cbind(imp.scanone[,1:2], apply(attr(imp.scanone, "lod.imp"), 1, FUN))
  attributes(output)=attributes(imp.scanone)
  attr(output, "FUN")=as.character(substitute(FUN))
  return(output)
} # update.imp 
