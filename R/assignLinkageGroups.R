############################################################################################################
#
# assignLinkageGroups.R
#
# Copyright (c) 2011, Danny Arends
#
# Modified by Konrad Zych
# 
# first written October 2011
# last modified December 2011
# last modified in version: 0.9.1-0
# in current version: active, in main workflow
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
# Contains: assignLinkageGroups, regorganizeMarkersWithin
#
############################################################################################################

############################################################################################################
#									*** assignLinkageGroups ***
#
# DESCRIPTION:
#	Assign linkage groups based on a user supplied known number of chromosomes
# we can use the genetic map of a cross, or the rf matrix
# 
# PARAMETERS:
#	cross - an object of class cross
#	n.chr - number of chromosomes expected
#	use - what kind of data should be used:
#		geno - genotypes of cross (cross$geno)
#		rf - recombination fractions (cross$rf)
#
# OUTPUT:
#	 an object of class cross
#
############################################################################################################
assignLinkageGroups <- function(cross, n.chr, use=c("geno","rf"), ...){
  inListCheck.internal(use,"use",c("rf","geno"))
  if(use=="geno") clustering <- kmeans(t(pull.geno(cross)), n.chr, ...)
  if(use=="rf") clustering <- kmeans(est.rf(cross)$rf, n.chr, ...)
  reorganizeMarkersWithin(cross, clustering$cluster)
}

############################################################################################################
#									*** regorganizeMarkersWithin ***
#
# DESCRIPTION:
#	Function to quickly rearange all the markers in a cross based on any ordering
# 
# PARAMETERS:
#	cross - an object of class cross
#	ordering - vector specifying new ordering of the markers
#
# OUTPUT:
#	 an object of class cross
#
############################################################################################################
reorganizeMarkersWithin <- function(cross, ordering){
  cross <- clean(cross)
  n.markers <- nmar(cross)
  chrtype <- rep(sapply(cross$geno, class), n.markers)
  crosstype <- class(cross)[1]
  g <- pull.geno(cross)
  newChromosomes <- sort(unique(ordering))
  n.newChromosomes <- length(newChromosomes)
  cross$geno <- vector("list", n.newChromosomes)
  names(cross$geno) <- newChromosomes
  for (i in 1:n.newChromosomes) {
      selectedMarkers <- which(ordering == newChromosomes[i])
      cross$geno[[i]]$data <- g[, selectedMarkers, drop = FALSE]
      cross$geno[[i]]$map <- seq(0, by = 10, length = length(selectedMarkers))
      if (crosstype == "4way") {
          cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map, cross$geno[[i]]$map)
          colnames(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
      }
      else names(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
      thechrtype <- unique(chrtype[which(ordering == newChromosomes[i])])
      if (length(thechrtype) > 1) 
          warning("Problem with linkage group ", i, ": A or X?\n", 
            paste(thechrtype, collapse = " "))
      else class(cross$geno[[i]]) <- thechrtype
  }
  cross <- est.rf(cross)
  return(cross)
}
