#
# power.plot.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified January, 2013
# first written January, 2013
# Contains: power.plot.R
#

# cross.saturate
#
# DESCRIPTION:
#  Saturate an existing genetic map by adding markers derived from expression
# OUTPUT:
#  An object of class cross
#
power.plot <- function(cross1,cross2,qtl.thr=5,n.pheno=500,...){
  if(missing(cross1)){
    stop("No object of class cross.")
  }
  if(missing(cross2)){
    stop("No object of class cross.")
  }
  if (!any(class(cross1) == "cross")){
    stop("Input should have class \"cross\".")
  }
  if (!any(class(cross2) == "cross")){
    stop("Input should have class \"cross\".")
  }
  if(!is.numeric(qtl.thr) || qtl.thr<0){
    stop("qtl.thr should be a numeric value bigger than 0.")
  }
  if(nphe(cross1)!=nphe(cross2)){
     stop("Difference in phenotype number between crosses.")
  }
  if(!is.numeric(n.pheno) || n.pheno<0 || n.pheno>nphe(cross1)){
    stop("qtl.thr should be a numeric value between 0 and nphe(cross1).")
  }
  markers <- sample(1:nphe(cross1),n.pheno)
  res1 <- scanone(cross1,pheno.col=markers,...)
  res2 <- scanone(cross2,pheno.col=markers,...)
  maxes1 <- apply(res1[,3:ncol(resO)],2,max)
  maxes2 <- apply(res2[,3:ncol(resS)],2,max)
  plot(maxes1,maxes2,pch=20,xlab="LOD scores in the original cross",ylab="LOD scores in the saturated cross")
  abline(0,1,col="red")
  abline(v=qtl.thr,lty=2,col="grey")
  abline(h=qtl.thr,lty=2,col="grey")
  significant<-which(maxes1>qtl.thr)
  cat(sum(maxes1[significant]<maxes2[significant]),"which is:",sum(maxes1[significant]<maxes2[significant])/length(significant),"% of significant qtls show increase in power\n")
  invisible(sum(maxes1[significant]<maxes2[significant])/length(significant))
}
