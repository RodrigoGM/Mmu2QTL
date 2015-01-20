#source('C:/Documents and Settings/RJGularte/My Documents/LAB FILES/Ricardo/qtl.r')
## Replicated CIM v.2 
##      changes made with help of Eva Chan <eva.chan@esiro.au>
## r.cim fucntion will replicate the cim a set number of X times, and return the average lod score for the replicates
## times= number of times to be replicated.  default is 10
## cross:  object cross
## save: if save=TRUE, a data frame containing all the lod scores for the replicates will be saved in a data frame
##       called rcim.results
r.cim<- function(cross, pheno.col=1, times=10, n.marcovar=3, window=2, method="hk", save=FALSE, ...) {
require("qtl") || stop("qtl package not available.")

    getrownames<- rownames(cim(cross, pheno.col=pheno.col, method=method, n.marcovar=1,...))
    output<- replicate(times, cim(cross, pheno.col=pheno.col, method=method, n.marcovar=n.marcovar, ...))
#output<- replicate(times, cim(cross, pheno.col=pheno.col, method=method, n.marcovar=n.marcovar, ...))

  chr<-output[,1]$chr
  pos<-output[,1]$pos
  get.lods<-sapply(1:times, function(o) output[,o]$lod)
  lod<- apply(get.lods, 1, median)  
#  lod<<- apply(matrix(lapply(get.lods,cbind)), 1, mean)

if(save==TRUE) {  
  rcim.results<-data.frame(chr, pos, get.lods)
  rownames(rcim.results)= getrownames
  colnames(rcim.results)=c("chr", "pos", paste("Sim.", seq(times), sep=""))
  rcim.results<<-rcim.results
  }
  
  results=data.frame(chr, pos, lod)
  rownames(results)=getrownames                       

class(results)=c("scanone","data.frame")
attr(results, "method") = method
attr(results, "map") = pull.map(cross)

return(results)
} #r.cim


## mscan.rcim will perform a multiple trait scan on an r.cim function.
mscan.rcim<- function(cross, phenotypes, times=2, n.marcovar=3, window=2, method="hk", ...){
require("qtl")|| stop("qtl package not available.")

nphenos=length(phenotypes)
output=list(0)

for (i in 1:nphenos) {
   p=phenotypes[i]
   output[[i]]=r.cim(cross=cross, pheno.col=find.pheno(cross, p), times=times, n.marcovar=n.marcovar, window=window,
                 method=method, ...)
   names(output)[length(output)]=p
   }
class(output)="mscanone"
attr(output, "method")= method
attr(output, "map")= pull.map(cross)

return(output)
}  ##mscan.rcim


## mcim will return a scanone output with multiple lod columns corresponding to the phenotypes provided

mcim<- function(cross, phenotypes=c(1,2), n.marcovar=3, window=2, method='hk', 
        imp.method='imp', error.prob=0.0001, ...) {
require(qtl) || stop("qtl package not available.")

if (length(phenotypes)<= 1) {stop("you only have one phenotype, please use 'cim'")}

if (length(phenotypes)>= 2) {
        
      cim..<- as.data.frame(lapply(phenotypes, function(x) cim(cross, pheno.col=x, n.marcovar=n.marcovar, 
            window=window, method=method, error.prob=error.prob, ...)))
    
      cim..<- cim..[,c('chr','pos', 'lod', paste('lod', seq(length(phenotypes)-1), sep="."))]
      colnames(cim..)<- c('chr', 'pos', phenotypes)
      }
      
class(cim..)<- c('scanone', 'data.frame')
attr(cim.., "method")= method
attr(cim.., "map")= pull.map(cross)

return(cim..)

} ## mcim


## mr.cim  multimple phenotypes r.cim
mr.cim<- function(cross, phenotypes=c(1,2), times=100, n.marcovar=3, window=2, method='hk',
            imp.method='imp', error.prob=0.0001, save=FALSE, ...) {
require(qtl) || stop('qtl package not avialable.')

if (save!= FALSE) stop("Fool, I have not yet figured out how to \n save all the replications in more than one phenotype")

if (length(phenotypes)<= 1) {stop('you only have one phenotype, plese use "r.cim"')}

if (length(phenotypes)>= 2) {

      rcim..<- as.data.frame(lapply(phenotypes, function(x) r.cim(cross, pheno.col=x, n.marcovar=n.marcovar, 
            times=times, window=window, method=method, error.prob=error.prob, ...)))
    
      rcim..<- rcim..[,c('chr','pos', 'lod', paste('lod', seq(length(phenotypes)-1), sep="."))]
      colnames(rcim..)<- c('chr', 'pos', phenotypes)
        
      }

class(rcim..)<- c('scanone', 'data.frame')
attr(rcim.., "method")= method
attr(rcim.., "map")= pull.map(cross)

return(rcim..)


    } ## End mr.cim






###################using lapply
#> mscan.rcim<- function(cross, phenotypes, times=2, n.marcovar=3, window=2, method="hk", ...){
#+ require("qtl")|| stop("qtl package not available.")
#+
#+ output=list(0)
#+ output[[phenotype]]=lapply(phenotypes, function(p){ r.cim(cross=cross, pheno.col=find.pheno(cross, p), times=times, method=method, ...)
#+                           names(output)[length(output)]=p
#+                           })
#+ class(output)="mscanone"
#+ attr(output, "method")= method
#+ attr(output, "map")= pull.map(cross)
#+
#+ return(output)
#+ }  ##mscan.rcim



