############################################################################################################
#
# assignLinkageGroups.R
#
# Copyright (c) 2011, Danny Arends
#
# 
# first written Oct 2011
# last modified Oct 2011
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
#	assigns linkage groups based on a user supplied known number of chromosomes, we can use the genetic map 
#	of a cross, or the rf matrix
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
  if(use=="geno") cl <- kmeans(t(pull.geno(cross)), n.chr, ...)
  if(use=="rf")cl <- kmeans(est.rf(cross)$rf, n.chr, ...)
  reorganizeMarkersWithin(cross, cl$cluster)
}

############################################################################################################
#									*** regorganizeMarkersWithin ***
#
# DESCRIPTION:
#	function that is missing from R/qtl to quickly rearange all the markers in a cross based on any ordering
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
  ingrp <- ordering
  tab <- sort(table(ingrp), decreasing = TRUE)
  u <- names(tab)
  revgrp <- ingrp
  for (i in seq(along = u)) revgrp[ingrp == u[i]] <- i
  cross <- clean(cross)
  n.mar <- nmar(cross)
  tot.mar <- totmar(cross)
  chrtype <- rep(sapply(cross$geno, class), n.mar)
  crosstype <- class(cross)[1]
  g <- pull.geno(cross)
  cross$geno <- vector("list", max(revgrp))
  names(cross$geno) <- 1:max(revgrp)
  for (i in 1:max(revgrp)) {
      cross$geno[[i]]$data <- g[, revgrp == i, drop = FALSE]
      cross$geno[[i]]$map <- seq(0, by = 10, length = tab[i])
      if (crosstype == "4way") {
          cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map, cross$geno[[i]]$map)
          colnames(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
      }
      else names(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
      thechrtype <- unique(chrtype[revgrp == i])
      if (length(thechrtype) > 1) 
          warning("Problem with linkage group ", i, ": A or X?\n", 
            paste(thechrtype, collapse = " "))
      else class(cross$geno[[i]]) <- thechrtype
  }
  cross <- est.rf(cross)
  return(cross)
}
