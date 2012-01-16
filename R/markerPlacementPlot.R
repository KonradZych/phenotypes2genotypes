############################################################################################################
#
# markerPlacementPlot.R
#
# Copyright (c) 2011, Danny Arends
#
# Modified by Konrad Zych
# 
# first written November 2011
# last modified January 2012
# last modified in version: 1.0.0
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
#
#     This program is distributed in the hope that it will be useful
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
#
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: markerPlacementPlot
#
############################################################################################################

markerPlacementPlot <- function(population, placeUsing=c("qtl","correlation"),cross){
  if(missing(population)) stop("Please provide a population object\n")
  if(is.null(population$offspring$genotypes$real)){
    stop("No original genotypes in population$offspring$genotypes$real, load them in using intoPopulation\n")
  }
  check.population(population)
  if(placeUsing=="correlation"){
    if(missing(cross)){
      cross <- genotypesToCross.internal(population, "simulated")
    }
    cormatrix <- map2mapCorrelationMatrix(cross,population)
    s <- NULL
    p <- seq(1,5,0.1)
    for(x in p){
      maximums <- apply(abs(cormatrix),2,max)
      means <- apply(abs(cormatrix),2,mean)
      sds <- apply(abs(cormatrix),2,sd)
      selected <- which(maximums > (means+x*sds))
      s <- c(s,length(selected))
    }
  plot(p,s,type='o',main="Number of markers placed",xlab="corThreshold",ylab="# of markers")
  }else{
    if(is.null(population$offspring$genotypes$qtl)) stop("No qtl data in population$offspring$genotypes$qtl, run scanQTLs function first.")
    singleqtl <- NULL
    noqtl <- NULL
    multipleqtl <- NULL
    thrRange <- seq(0,20,0.1)
    for(tr in thrRange){
      peaks <- getpeaks.internal(population$offspring$genotypes$qtl,tr)
      rownames(peaks) <- rownames(population$offspring$genotypes$qtl)
      colnames(peaks) <- colnames(population$offspring$genotypes$qtl)
      nqtl <- apply(peaks,1,function(row){sum(row==2)})
      singleqtl <- c(singleqtl,sum(nqtl==1))
      noqtl <- c(noqtl,sum(nqtl==0))
      multipleqtl <- c(multipleqtl,sum(nqtl>1))
    }
    plot(thrRange,noqtl,type='o',col="red",main="Number of markers placed",xlab="threshold",ylab="# of markers",ylim=c(0,nrow(population$offspring$genotypes$qtl)))
    points(thrRange,singleqtl,type='o',col="green")
    points(thrRange,multipleqtl,type='o',col="blue")
    legend(x="topright",legend=c("no peak","single peak","multiple peaks"),col=c("red","green","blue"),cex=0.8,pch=21,lwd=2,bg="white")
  }
}
