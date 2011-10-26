############################################################################################################
#
# createNewMap.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified October 2011
# last modified in version: 0.9.1
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
# Contains: orderChromosomes, majorityRule.internal, mergeChromosomes.internal,
#			switchChromosomes.internal, removeChromosomes.internal, removeChromosomesSub.internal,
#
############################################################################################################

############################################################################################################
#									*** createNewMap ***
#
# DESCRIPTION:
# 	function to create new map and save it in cross object
# 
# PARAMETERS:
# 	population - object of class population
# 	orde - object of class population
# 	n.chr - expected number of linkage groups
# 	use - expected number of linkage groups
#	verbose - be verbose
#
# OUTPUT:
#	an object of class cross
#
#
############################################################################################################
createNewMap <- function(population,  n.chr,  use=c("geno","rf"), verbose=FALSE, debugMode=0){
	if(missing(n.chr)) stop("n.chr in an obligatory parameter")
	if(missing(population)) stop("no population object provided")
	use <- defaultCheck.internal(use, "use", 2, "geno")
	check.population(population)
	if(is.null(population$offspring$genotypes$simulated)){
		stop("no simulated genotypes in population object, first use toGenotypes!\n")
	}
	#*******SAVING CROSS OBJECT*******
	s1 <- proc.time()
  aa <- tempfile()
	sink(aa)
	cross <- genotypesToCross.internal(population,"simulated",verbose=verbose,debugMode=debugMode)
  sink()
	file.remove(aa)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("saving data into cross object done in:",(e1-s1)[3],"seconds.\n")
	
	#####
	######
	########## VERY DIRTY HACK
	cross <- convert2riself(cross)
	
	#*******CREATING NEW MAP*******
	s1 <- proc.time()
	cross <- assignLinkageGroups(cross,n.chr)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("New map created in:",(e1-s1)[3],"seconds.\n")
	
	invisible(cross)
}