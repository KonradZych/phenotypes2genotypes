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

markerPlacementPlot <- function(population){
  mycross <- genotypesToCross.internal(population, "simulated")
  cormatrix <- map2mapCorrelationMatrix(mycross,population)
  s <- NULL
  p <- seq(0.1,1,0.01)
  for(x in p){
    s <- c(s,length(which(apply(abs(cormatrix),2,max)>x)))
  }
  plot(p,s,type='o',main="Number of markers placed",xlab="corThreshold",ylab="# of markers")
}
