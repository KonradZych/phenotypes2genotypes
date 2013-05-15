#
# power.plot.R
#
# Copyright (c) 2010-2013 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2013
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
power.plot <- function(cross1,cross2,qtlThr=5,nPheno=500,verbose=FALSE,...){
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
  if(!is.numeric(qtlThr) || qtlThr<0){
    stop("qtlThr should be a numeric value bigger than 0.")
  }
  if(nphe(cross1)!=nphe(cross2)){
     stop("Difference in phenotype number between crosses.")
  }
  if(!is.numeric(nPheno) || nPheno<2 || nPheno>nphe(cross1)){
    stop("nPheno should be a numeric value between 2 and nphe(cross1).")
  }
  if(is.null(cross1$geno[[1]]$prob)){
    cross1 <- calc.genoprob(cross1)
  }
  if(is.null(cross2$geno[[1]]$prob)){
    cross2 <- calc.genoprob(cross2)
  }
  markers  <- sample(1:nphe(cross1),nPheno)
  qtlThr2 <- qtlThr
  res1     <- scanone(cross1,pheno.col=markers,...)
  res2     <- scanone(cross2,pheno.col=markers,...)
  maxes1   <- apply(res1[,3:ncol(res1)],2,max)
  maxes2   <- apply(res2[,3:ncol(res2)],2,max)
  if(verbose){
    vals           <- cbind(markers,maxes1,maxes2)
    colnames(vals) <- c("#phenotype","max LOD in cross1","max LOD in cross2")
    print(vals)
  }
  colorNr   <- as.numeric(maxes1>maxes2)+1
  colorCols <- c("green","red")
  plot(maxes1,maxes2,pch=20,xlab="LOD scores in the original cross",ylab="LOD scores in the saturated cross",main="Figure 3 - QTL detection power.",col=colorCols[colorNr])
  abline(qtlThr2-qtlThr,1,col="red") #diagonal line
  abline(v=qtlThr,lty=2,col="grey")   #dotted line showing threshold used for x axis
  abline(h=qtlThr2,lty=2,col="grey")  #dotted line showing threshold used for y axis
  significant  <- which(maxes1>qtlThr)  #nr of significant QTLs in cross 1
  significant2 <- which(maxes2>qtlThr2) #nr of significant QTLs in cross 2
  stronger     <- sum(maxes1[significant]<(maxes2[significant])) #nr of QTLs showubg increase in power
  cat(length(significant),"significant qtls\n")
  cat(stronger,"which is:",stronger/length(significant),"% of significant qtls show increase in power\n")
  cat(length(significant2)-length(significant),"new QTLs found\n") # nr of gained QTLs
  invisible(list(res1,res2))
}

