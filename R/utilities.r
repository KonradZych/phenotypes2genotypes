############################################################################################################
#
# utilities.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified June 2011
# last modified in version: 0.8.1
# in current version: active, not in main workflow
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
#
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
#
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: print.population, removeIndividuals, doCleanUp.internal
#
############################################################################################################

############################################################################################################
#									*** print.population ***
#
# DESCRIPTION:
#	overwrites the print function for objects of class "population"
# 
# PARAMETERS:
# 	x - object of class population
# 	... - passed to cats
# 
# OUTPUT:
#	none
#
############################################################################################################
print.population <- function(x,...){
	cat("This is object of class \"population\"\n  It is too complex to print, so we provide just this summary.\n")
	if(!(is.null(x$offspring))){
    cat("Offspring:\n",...)
		if(!(is.null(x$offspring$phenotypes))){
			cat("\tPhenotypes:",ncol(x$offspring$phenotypes),"\n",...)
      cat("\tMarkers:",nrow(x$offspring$phenotypes),"\n",...)
		}else{
			stop("No phenotype data for offspring, this is not a valid population object\n")
		}
		if(!(is.null(x$offspring$genotypes$read))){
			cat("\tGenotypes:",ncol(x$offspring$genotypes),"\n",...)
		}else{
			cat("\tGenotypes: None\n",...)
		}
		if(!(is.null(x$maps$genetic))){
			cat("\tGenetic map:",nrow(x$maps$genetic),"markers, ",length(table(x$maps$genetic[,1]))," chromosomes\n",...)
		}else{
			cat("\ttGenetic map: None\n")
		}
		if(!(is.null(x$maps$physical))){
			cat("\tPhysical map:",nrow(x$maps$physical),"markers, ",length(table(x$maps$physical[,1]))," chromosomes\n",...)
		}else{
			cat("\tPhysical ap: None\n")
		}    
	}else{
		stop("No phenotype data for offspring, this is not a valid population object\n")
	}
	
	if(!(is.null(x$founders))){
    cat("Founders:\n",...)
		if(!(is.null(x$founders$phenotypes))){
			cat("\tPhenotypes:",ncol(x$founders$phenotypes),"\n",...)
      cat("\tMarkers:",nrow(x$founders$phenotypes),"\n",...)
		}else{
			stop("No phenotype data for founders, this is not a valid population object\n")
		}
		if(!(is.null(x$founders$RP))){
			cat("\tDifferential expression: Detected\n",...)
		}else{
			cat("\tDifferential expression: Not Detected (please: use functionname) \n",...)
		}
		if(!(is.null(x$founders$groups))){
			cat("\tFounder groups:",x$founders$groups,"\n",...)
		}else{
			stop("No information about founders groups\n",...)
		}
	}else{
		stop("No phenotype data for founders, this is not a valid population object\n")
	}
}

############################################################################################################
#									*** removeIndividuals ***
#
# DESCRIPTION:
#	Function to remove individual(s) from population object. 
# 
# PARAMETERS:
# 	population - object of class population
# 	individuals - individuals to be romved specified by their names
#
# OUTPUT:
#	object of class population
#
############################################################################################################
removeIndividuals <- function(population,individuals,verbose=FALSE){
	for(ind in individuals){
		if(ind%in%colnames(population$offspring$genotypes$real)){
			population$offspring$genotypes$real <- population$offspring$genotypes$real[,-which(colnames(population$offspring$genotypes$real)==ind)]
			if(verbose)cat("Removed",ind,"from population$offspring$genotypes$real\n")
		}
		if(ind%in%colnames(population$offspring$phenotypes)){
			population$offspring$phenotypes <- population$offspring$phenotypes[,-which(colnames(population$offspring$phenotypes)==ind)]
			if(verbose)cat("Removed",ind,"from population$offspring$phenotypes\n")
		}
		if(ind%in%colnames(population$founders$phenotypes)){
			population$founders$phenotypes <- population$founders$phenotypes[,-which(colnames(population$founders$phenotypes)==ind)]
			if(verbose)cat("Removed",ind,"from population$founders$phenotypes\n")
		}
	}
	invisible(population)
}


############################################################################################################
#									*** doCleanUp.internal ***
#
# DESCRIPTION:
#	better garbage collection 
# 
# PARAMETERS:
#	verbose - be verbose
#
# OUTPUT:
#	none
#
############################################################################################################
doCleanUp.internal <- function(verbose=FALSE){
	before <- gc()[2,3]
	bf <- before
	after <- gc()[2,3]
	while(before!=after){
		before <- after
		after <- gc()[2,3]
	}
	if(verbose) cat("Cleaned up memory from:",bf,"to:",after,"\n")
}

