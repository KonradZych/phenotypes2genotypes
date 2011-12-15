############################################################################################################
#
# markerPlacementPlot.R
#
# Copyright (c) 2011, Danny Arends
#
# Modified by Konrad Zych
# 
# first written November 2011
# last modified December 2011
# last modified in version: 0.9.1
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

markerPlacementPlot <- function(population, placeUsing=c("qtl","correlation"), cross){
  if(missing(cross)){
    cross <- genotypesToCross.internal(population, "simulated")
  }
  if(placeUsing=="correlation"){
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
    genotypes <- population$offspring$genotypes$real
    markers <- markernames(cross)
    phenotypes <- pull.geno(cross)[,markers]
    results <- NULL
    rnames <- NULL
    cnt <- 1
    for(x in 1:ncol(phenotypes)){
      results <- rbind(results,apply(genotypes,1, 
        function(geno){
          linmod <- lm(phenotypes[,x] ~ geno)
          -log10(anova(linmod)[[5]][1])
        }
      ))
      rnames <- c(rnames,colnames(phenotypes)[x])
      cnt <- cnt +1
    }
  rownames(results) <- rnames
  singleqtl <- NULL
  noqtl <- NULL
  multipleqtl <- NULL
  for(tr in seq(0,20,0.5)){
    nqtls <- apply(results,1,function(x){sum(x>tr)})
    singleqtl <- c(singleqtl,sum(nqtls==1))
    noqtl <- c(noqtl,sum(nqtls==0))
    multipleqtl <- c(multipleqtl,sum(nqtls>1))
  }
  plot(noqtl,type='l',col="red",main="Number of markers placed",xlab="threshold",ylab="# of markers")
  points(singleqtl,type='l',col="green")
  points(multipleqtl,type='l',col="blue")
  }
}
