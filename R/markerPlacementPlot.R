############################################################################################################
#
# markerPlacementPlot.R
#
# Copyright (c) 2011, Danny Arends
#
# 
# first written Nov 2011
# last modified Nov 2011
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

markerPlacementPlot <- function(population, cross){
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
}
