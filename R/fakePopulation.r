############################################################################################################
#
# fakePopulation.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified September 2011
# last modified in version: 0.9.0
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
# Contains: fakePopulation,
#             fakePheno.internal, fakeFounders.internal, fakePhysicalMap.internal, convertMap.internal, 
#             simBC.internal, fakeMixUps
#
############################################################################################################


############################################################################################################
#                                          *** fakePopulation ***
#
# DESCRIPTION:
#	simulating object of class population 
# 
# PARAMETERS:
#	n.founders - number of founders to be simulated
#	n.offspring - number of offspring individuals to be simulated
#	n.markers - number of markers to be simulated
#	n.chromosomes - number of chromosomes to be simulated

#	type - type of the cross 
#		- "riself" - RIL by selfing
#		- "f2" - f2 cross
#		- "bc" - back cross
#		- "risib" - RIL by sibling mating
#	n.mixups - number of mixups to be simulated in data
#	... - parameters send to sim.cross function
#
# OUTPUT:
#	object of class population
#
############################################################################################################
fakePopulation <- function(n.founders = 4, n.offspring = 100, n.markers=100,n.chromosomes=10, type = c("riself", "f2", "bc", "risib"), n.mixups=0, verbose=FALSE,...){
  type <- match.arg(type)
  if(!(is.numeric(n.founders))) stop("n.founders should be numeric\n")
  if(!(is.numeric(n.offspring))) stop("n.offspring should be numeric\n")
  if(!(is.numeric(n.markers))) stop("n.markers should be numeric\n")
  if(!(is.numeric(n.chromosomes))) stop("n.chromosomes should be numeric\n")
  if(n.founders<4){
    warning("too small n.founders, changing to 4\n")
    n.founders <- 4
  }
  
	if(!(n.founders%%2==0)){
    warning("n.founders should be even, changing to",n.founders+1,"\n")
    n.founders <- n.founders+1
   }
	if(length(type)>1) type <- type[1]
	if(n.offspring<10){
		warning("too small n.offspring, changing to 10\n")
		n.offspring <- 10
	}
	if(n.markers<10){
		warning("too small n.markers, changing to 10\n")
		n.markers <- 10
	}
	if(n.markers<n.chromosomes){
		warning("n.markers cannot be smaller than n.chromosomes, changing n.markers to 10*n.chromosomes\n")
		n.markers <- 10*n.chromosomes
	}

	### fuction itself
	map <- sim.map(rep(100,n.chromosomes),n.mar=(n.markers/n.chromosomes), include.x=FALSE,)
	fake <- sim.cross(map,type=type, n.ind=n.offspring, ...)
	geno <- t(pull.geno(fake))
	if(type=="bc"){
		geno <- t(apply(geno,1,simBC.internal))
	}
	map <- convertMap.internal(map)
	physicalMap <- fakePhysicalMap.internal(map)
	colnames(geno) <- paste("RIL",1:ncol(geno),sep="_")
	pheno <- t(apply(geno,1,fakePheno.internal))
	rownames(pheno) <- rownames(geno)
	colnames(pheno) <- colnames(geno)
  if(n.mixups>0){
    pheno <- fakeMixUps.internal(pheno,n.mixups)
  }
	founders <- t(apply(pheno,1,fakeFounders.internal,n.founders))
	rownames(founders) <- rownames(geno)
	colnames(founders) <- 1:n.founders
	colnames(founders)[1:(n.founders/2)] <- paste("Founder",1,1:(n.founders/2),sep="_")
	colnames(founders)[(n.founders/2+1):n.founders] <- paste("Founder",2,(n.founders/2+1):n.founders,sep="_")
	geno[which(geno==2)] <- 0
	foundersGroups <- c(rep(0,(n.founders/2)),rep(1,(n.founders/2)))
	population <- createPopulation(pheno, founders, foundersGroups, geno, map, physicalMap,verbose=verbose)
	check.population(population)
	invisible(population)
}

############################################################################################################
#                                  *** fakePheno.internal ***
#
# DESCRIPTION:
#	subfunction of fakePopulation - simulating phenotype data using genotype data simulated by sim.cross
# 
# PARAMETERS:
#	genoRow - row of offspring genotype matrix
#
# OUTPUT:
#	row of offspring phenotype matrix
#
############################################################################################################
fakePheno.internal <- function(genoRow,maxScale=10,maxError=3){
	scalingF <- runif(1,1,maxScale)
	errorF <- runif(length(genoRow),0,maxError)
	genoRow <- (genoRow*scalingF) + errorF
	invisible(genoRow)
}

############################################################################################################
#                                  *** fakeFounders.internal ***
#
# DESCRIPTION:
#	subfunction of fakePopulation - simulating founders phenotype data using offspring phenotype data
# 
# PARAMETERS:
#	phenoRow - row of offspring phenotype matrix
#
# OUTPUT:
#	row of parental phenotype matrix
#
############################################################################################################
fakeFounders.internal <- function(phenoRow,n.founders){
	errorF <- runif(n.founders,0,2)
	up <- runif(1,-1,1)
	diffExprRate <- runif(1,0.1,1)
	cur_mean <- mean(phenoRow)
	if(up>=0){
		foundersRow <- c(rep((cur_mean-diffExprRate*cur_mean),(n.founders/2)),rep(cur_mean+diffExprRate*cur_mean,(n.founders/2))) + errorF
	}else{
		foundersRow <- c(rep((cur_mean+diffExprRate*cur_mean),(n.founders/2)),rep(cur_mean-diffExprRate*cur_mean,(n.founders/2))) - errorF
	}
	invisible(foundersRow)
}

############################################################################################################
#                                  *** fakePhysicalMap.internal ***
#
# DESCRIPTION:
#	simulating physical map using genetic one
# 
# PARAMETERS:
#	map - map as used in population object (matrix with two cols 1-chr nr, 2- position)
#
# OUTPUT:
#	map of the same type
#
############################################################################################################
fakePhysicalMap.internal <- function(map){
	for(i in 1:nrow(map)){
		errorF <- runif(1,0,100)
		if(errorF>90){
			newChrom <- round(runif(1,1,(max(map[,1]))))
			map[i,1] <- newChrom
		}
	}
	invisible(map)
}

############################################################################################################
#                                  *** convertMap.internal ***
#
# DESCRIPTION:
#	convert rqtl type map into population type one
# 
# PARAMETERS:
#	map - map of cross type (list with names - names of chromosomes, elements o the list - markers and their
#		  positions)
#
# OUTPUT:
#	map of population class type - rownames - names of the markers, first column - numbers of chromosomes,
#		second - position of the marker on chromosome
#
############################################################################################################
convertMap.internal <- function(map){
	map_ <- NULL
	for(i in 1:length(map)){
		cur_chr <- cbind(rep(i,length(map[[i]])),map[[i]])
		map_ <- rbind(map_,cur_chr)
	}
	invisible(map_)
}

############################################################################################################
#                                  *** simBC.internal ***
#
# DESCRIPTION:
#	convert rqtl type map into population type one
# 
# PARAMETERS:
#	map - map of cross type (list with names - names of chromosomes, elements o the list - markers and their
#		  positions)
#
# OUTPUT:
#	map of population class type - rownames - names of the markers, first column - numbers of chromosomes,
#		second - position of the marker on chromosome
#
############################################################################################################
simBC.internal <- function(genoRow){
	n.ind <- length(genoRow)
	toCount <- round(runif(1,1,2))
	toChange <- (3-toCount)
	numberOfToCount <- sum(genoRow==toCount)
	errorRate <- round(n.ind*(runif(1,0,5)/100))
	expected <- round(0.75*n.ind)
	cat(numberOfToCount,expected,errorRate,"\n")
	while(numberOfToCount<(expected-errorRate)){
		genoRow[round(runif(1,1,length(genoRow)))] <- toCount
		numberOfToCount <- sum(genoRow==toCount)
	}
	while(numberOfToCount>(expected+errorRate)){
		genoRow[round(runif(1,1,length(genoRow)))] <- toChange
		numberOfToCount <- sum(genoRow==toCount)
	}
	other <- sum(genoRow==toChange)
	cat(numberOfToCount,other ,"\n")
	invisible(genoRow)
}

############################################################################################################
#                                  *** mixUpPheno.internal ***
#
# DESCRIPTION:
#	convert rqtl type map into population type one
# 
# PARAMETERS:
#	map - map of cross type (list with names - names of chromosomes, elements o the list - markers and their
#		  positions)
#
# OUTPUT:
#	map of population class type - rownames - names of the markers, first column - numbers of chromosomes,
#		second - position of the marker on chromosome
#
############################################################################################################
fakeMixUps.internal<- function(pheno, n.mixups){
  for(i in 1:n.mixups){
    toMix <- sample(1:nrow(pheno),2)
    temp_ <- pheno[toMix[1],]
    pheno[toMix[1],] <- pheno[toMix[2],]
    pheno[toMix[2],]<- pheno[toMix[1],]
  }
  invisible(pheno)
}
