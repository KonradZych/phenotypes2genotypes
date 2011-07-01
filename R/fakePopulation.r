############################################################################################################
#
# utilities.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified July 2011
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
# Contains: fakePopulation, fakePheno.internal, fakeFounders.internal, convertMap.internal
#
############################################################################################################


############################################################################################################
#									*** fakePopulation ***
#
# DESCRIPTION:
#	simulating population object 
# 
# PARAMETERS:
#	nrFounders - number of founders to be simulated
#	... - parameters send to sim.cross, most important:
#		type - type of the cross - "riself" for two-way RIL
#		n.ind - number of offspring individuals
#
# OUTPUT:
#	object of class population
#
############################################################################################################
fakePopulation <- function(nrFounders = 4,...){
	if(nrFounders<4) nrFounders <- 4
	if(!(nrFounders%%2==0)) nrFounders <- nrFounders+1
	map <- sim.map()
	if(!missing(...)){
		dots <- list(...)
		if ("type" %in% names(dots)){
			type <- dots$type
			#dots <- dots[[-which(names(dots)=="type")]]
		}else{
			type <- "riself"
		}
		if ("n.ind" %in% names(dots)){
			n.ind <- dots$n.ind
			#dots <- dots[[-which(names(dots)=="n.ind")]]
		}else{
			n.ind <- 250
		}
		if ("model" %in% names(dots)){
			model <- dots$model
			#dots <- dots[[-which(names(dots)=="model")]]
		}else{
			model <- rbind(c(1,45,1,1),c(5,20,0.5,-0.5))
		}
	}else{
		type <- "riself"
		n.ind <- 250
		model <- rbind(c(1,45,1,1),c(5,20,0.5,-0.5))
	}
	fake <- sim.cross(map,type=type, n.ind=n.ind, model = model)
	geno <- t(pull.geno(fake))
	map <- convertMap.internal(map)
	colnames(geno) <- paste("RIL",1:ncol(geno),sep="_")
	pheno <- t(apply(geno,1,fakePheno.internal))
	rownames(pheno) <- rownames(geno)
	colnames(pheno) <- colnames(geno)
	founders <- t(apply(pheno,1,fakeFounders.internal,nrFounders))
	rownames(founders) <- rownames(geno)
	colnames(founders) <- 1:nrFounders
	colnames(founders)[1:(nrFounders/2)] <- paste("Founder",1,1:(nrFounders/2),sep="_")
	colnames(founders)[(nrFounders/2+1):nrFounders] <- paste("Founder",2,(nrFounders/2+1):nrFounders,sep="_")
	geno[which(geno==2)] <- 0
	foundersGroups <- c(rep(0,(nrFounders/2)),rep(1,(nrFounders/2)))
	population <- createPopulation(pheno, founders, foundersGroups, geno, map, map)
	is.population(population)
	invisible(population)
}

############################################################################################################
#									*** fakePheno.internal ***
#
# DESCRIPTION:
#	simulating phenotype data using genotype data simulatyed by sim.cross
# 
# PARAMETERS:
#	genoRow - row of offspring genotype matrix
#
# OUTPUT:
#	row of offspring phenotype matrix
#
############################################################################################################
fakePheno.internal <- function(genoRow){
	scalingF <- runif(1,1,10)
	errorF <- runif(length(genoRow),0,3)
	genoRow <- (genoRow*scalingF) + errorF
	invisible(genoRow)
}

############################################################################################################
#									*** fakeFounders.internal ***
#
# DESCRIPTION:
#	simulating founders phenotype data using offspring phenotype data
# 
# PARAMETERS:
#	phenoRow - row of offspring phenotype matrix
#
# OUTPUT:
#	row of parental phenotype matrix
#
############################################################################################################
fakeFounders.internal <- function(phenoRow,nrFounders){
	errorF <- runif(nrFounders,0,3)
	cur_mean <- mean(phenoRow)
	foundersRow <- c(rep((cur_mean-0.1*cur_mean),(nrFounders/2)),rep((cur_mean-0.1*cur_mean),(nrFounders/2))) + errorF
	invisible(foundersRow)
}

############################################################################################################
#									*** convertMap.internal ***
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