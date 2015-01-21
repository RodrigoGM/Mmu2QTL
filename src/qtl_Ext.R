# PUTATIVE FUNCTIONS
mscanone<- function(cross, phenotypes, addcov=NULL, addname=NULL, intname=NULL, intocov=NULL, method="hk", ...) {

  lapply(phenotypes, function(pheno) scanone(cross, pheno.col=find.pheno(cross,pheno), method=method,  addcovar=adcov, intcovar=intcov, ...)
}
