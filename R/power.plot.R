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
power.plot <- function(cross1,cross2,qtl.thr=5,n.pheno=500,verbose=FALSE,...){
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
  if(!is.numeric(n.pheno) || n.pheno<2 || n.pheno>nphe(cross1)){
    stop("n.pheno should be a numeric value between 2 and nphe(cross1).")
  }
  if(is.null(cross1$geno[[1]]$prob)){
    cross1 <- calc.genoprob(cross1)
  }
  if(is.null(cross2$geno[[1]]$prob)){
    cross2 <- calc.genoprob(cross2)
  }
  markers <- sample(1:nphe(cross1),n.pheno)
  qtl.thr2 <- qtl.thr #- log10(sum(nmar(cross1))/sum(nmar(cross2)))
  print(qtl.thr2)
  res1 <- scanone(cross1,pheno.col=markers,...)
  res2 <- scanone(cross2,pheno.col=markers,...)
  maxes1 <- apply(res1[,3:ncol(res1)],2,max)
  maxes2 <- apply(res2[,3:ncol(res2)],2,max)
  if(verbose){
    vals <- cbind(markers,maxes1,maxes2)
    colnames(vals) <- c("#phenotype","max LOD in cross1","max LOD in cross2")
    print(vals)
  }
  plot(maxes1,maxes2,pch=20,xlab="LOD scores in the original cross",ylab="LOD scores in the saturated cross",main="Figure 3 - QTL detection power.")
  abline(qtl.thr2-qtl.thr,1,col="red")
  abline(v=qtl.thr,lty=2,col="grey")
  abline(h=qtl.thr2,lty=2,col="grey")
  significant<-which(maxes1>qtl.thr)
  stronger <- sum(maxes1[significant]<(maxes2[significant]))#+log10(sum(nmar(cross1))/sum(nmar(cross2)))))
  cat(length(significant),"significant qtls\n")
  cat(stronger,"which is:",stronger/length(significant),"% of significant qtls show increase in power\n")
  cat(sum(maxes2>qtl.thr2)-sum(maxes1>qtl.thr),"new QTLs found\n")
  invisible(list(res1,res2))
}

