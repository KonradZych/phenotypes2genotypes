############################################################################################################
#
# utilities.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified November 2011
# last modified in version: 0.9.1
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
# Contains: print.population, remove.individuals, doCleanUp.internal
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
print.population <- function(x, ...){
	if(missing(x)) stop("Please, provide an object of class population.\n")
	if(!(any(names(x)=="offspring"))){  stop("This is not correct population object.\n") }
	if(!(any(names(x)=="founders"))){  stop("This is not correct population object.\n") }
	cat("This is object of class \"population\"\n  It is too complex to print, so we provide just this summary.\n\n")
	if(class(x)[2] == "riself"){
		cat("Population type: RILs by selfing\n\n")
	}else if(class(x)[2] == "risib"){
		cat("Population type: RILs by sibling mating\n\n")
	}else if(class(x)[2] == "bc"){
		cat("Population type: backcross\n\n")
	}else if(class(x)[2] == "f2"){
		cat("Population type: F2 intercross\n\n")
	}
	if(!(is.null(x$offspring))){
    cat("Offspring (",ncol(x$offspring$phenotypes),"):\n",sep="",...)
		if(!(is.null(x$offspring$phenotypes))){
      cat("\tPhenotypes:",nrow(x$offspring$phenotypes),"\n",...)
		}else{
			stop("No phenotype data for offspring, this is not a valid population object\n")
		}
		if(!(is.null(x$offspring$genotypes$real))){
			cat("\tOriginal genotypes:",ncol(x$offspring$genotypes$real),"individuals",nrow(x$offspring$genotypes$real),"markers\n",...)
		}else{
			cat("\tOriginal genotypes: None\n",...)
		}
    if(!(is.null(x$offspring$genotypes$simulated))){
			cat("\tSimulated genotypes:",ncol(x$offspring$genotypes$simulated),"individuals",nrow(x$offspring$genotypes$simulated),"markers\n",...)
		}else{
			cat("\tSimulated genotypes: None\n",...)
		}
		if(!(is.null(x$offspring$genotypes$qtl))){
			cat("\tQTL scan results detected.\n",...)
		}else{
			cat("\tQTL scan results not detected, run scan.qtls.\n",...)
		}
		if(!(is.null(x$maps$genetic))){
			cat("\tGenetic map:",nrow(x$maps$genetic),"markers, ",length(table(x$maps$genetic[,1]))," chromosomes\n",...)
		}else{
			cat("\tGenetic map: None\n")
		}
		if(!(is.null(x$maps$physical))){
			cat("\tPhysical map:",nrow(x$maps$physical),"markers, ",length(table(x$maps$physical[,1]))," chromosomes\n",...)
		}else{
			cat("\tPhysical map: None\n")
		}    
	}else{
		stop("No phenotype data for offspring, this is not a valid population object\n")
	}

	if(!(is.null(x$founders))){
    cat("Founders (",ncol(x$founders$phenotypes),"):\n",sep="",...)
		if(!(is.null(x$founders$phenotypes))){
      cat("\tPhenotypes:",nrow(x$founders$phenotypes),"\n",...)
		}else{
			stop("No phenotype data for founders, this is not a valid population object\n")
		}
		if(!(is.null(x$founders$RP))){
			cat("\tDifferential expression: Detected\n",...)
		}else{
			cat("\tDifferential expression: Not Detected (please: use find.diff.expressed) \n",...)
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
#									*** remove.individuals ***
#
# DESCRIPTION:
#	Function to remove individual(s) from population object. 
# 
# PARAMETERS:
# 	population - object of class population
# 	individuals - vector of individuals to be removed specified by name
#
# OUTPUT:
#	object of class population
#
############################################################################################################
remove.individuals <- function(population,individuals,verbose=FALSE){
	check.population(population)
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
		if(ind%in%colnames(population$offspring$genotypes$simulated)){
			population$offspring$genotypes$simulated <- population$offspring$genotypes$simulated[,-which(colnames(population$offspring$genotypes$simulated)==ind)]
			if(verbose)cat("Removed",ind,"from population$offspring$genotypes$simulated\n")
		}
	}
	invisible(population)
}


############################################################################################################
#									*** doCleanUp.internal ***
#
# DESCRIPTION:
#	Force garbage collection in R, untill no objects can be cleaned anymore
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

############################################################################################################
#									*** write.population ***
#
# DESCRIPTION:
#	Writing population object to specific population file
# 
# PARAMETERS:
#	population - object of class population
#	outputFile - name of the output file
#	verbose - be verbose
#
# OUTPUT:
#	none
#
############################################################################################################
write.population <- function(population,outputFile="population.txt",verbose=FALSE){
	check.population(population)
	firstLine <- vector(mode="numeric",length=6)
	firstLine[1] <- nrow(population$offspring$phenotypes)
	firstLine[2] <- nrow(population$founders$phenotypes)
	firstLine[3] <- length(population$founders$groups)
	if(!is.null(population$offspring$genotypes$real)){
		firstLine[4] <- nrow(population$offspring$genotypes$real)
	}else{
		firstLine[4] <- NA
	}
	if(!is.null(population$maps$genetic)){
		firstLine[5] <- nrow(population$maps$genetic)
	}else{
		firstLine[5] <- NA
	}
	if(!is.null(population$maps$physical)){
		firstLine[6] <- nrow(population$maps$physical)
	}else{
		firstLine[6] <- NA
	}
	cat(firstLine,"\n",file=outputFile,append=FALSE)

	write.table(population$offspring$phenotypes,file=outputFile,sep="\t",quote=FALSE,append=TRUE)
	if(verbose) cat("Offspring phenotype data written into",outputFile,"\n")

	write.table(population$founders$phenotypes,file=outputFile,sep="\t",quote=FALSE,append=TRUE)
	if(verbose) cat("Founders phenotype data written into",outputFile,"\n")

	cat(population$founders$groups,"\n",file=outputFile,append=TRUE)
	if(verbose) cat("Information about founders groups data written into",outputFile,"\n")

	if(!is.null(population$offspring$genotypes$real)){
		write.table(population$offspring$genotypes$real,file=outputFile,sep="\t",quote=FALSE,append=TRUE)
		if(verbose) cat("Offspring genotype data written into",outputFile,"\n")
	}else{
		if(verbose) cat("Offspring genotype data not found.\n")
	}

	if(!is.null(population$maps$genetic)){
		write.table(population$maps$genetic,file=outputFile,sep="\t",quote=FALSE,append=TRUE,col.names=FALSE)
		if(verbose) cat("Genetic map written into",outputFile,"\n")
	}else{
		if(verbose) cat("Genetic map not found.\n")
	}

	if(!is.null(population$maps$physical)){
		write.table(population$maps$physical,file=outputFile,sep="\t",quote=FALSE,append=TRUE,col.names=FALSE)
		if(verbose) cat("Physical map written into",outputFile,"\n")
	}else{
		if(verbose) cat("Physical map not found.\n")
	}
}

############################################################################################################
#                                          ** assignChrToMarkers***
#
# DESCRIPTION:
#   Creating ordering vector from chromosomes assignment vector
# OUTPUT:
#  Vector for each of the markers specifying into which chromosome it should be moved
#
############################################################################################################
assignChrToMarkers <- function(assignment,cross){
    ordering <- vector(sum(nmar(cross)),mode="numeric")
    names(ordering) <- markernames(cross)
    for(i in 1:length(assignment)){
      oldChrom <- as.numeric(names(assignment)[i])
      newChrom <- assignment[i]
      markersFromOldChrom <- colnames(cross$geno[[oldChrom]]$data)
      ordering[markersFromOldChrom] <- rep(newChrom,length(markersFromOldChrom))
    }
    invisible(ordering)
}

############################################################################################################
#                                          ** pull.geno.from.cross***
#
# DESCRIPTION:
#   Pulling genotypes with a map from cross and putting into population object.
# OUTPUT:
#  An object of class population
#
############################################################################################################
set.geno.from.cross <- function(cross,population,map=c("genetic","physical")){
  map <- checkParameters.internal(map,c("genetic","physical"),"map")
  if(missing(population)) stop("Please provide a population object\n")
  if(missing(cross)) stop("Please provide a cross object\n")
  check.population(population)
  population$offspring$genotypes$real <- t(pull.geno(cross))
  if(map=="genetic"){
    population$maps$genetic <- convertMap.internal(pull.map(cross))
  }else{
    population$maps$physical <- convertMap.internal(pull.map(cross))
  }
  invisible(population)
}
