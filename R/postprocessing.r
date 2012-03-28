############################################################################################################
#
# postprocessing.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified November 2011
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
#									*** mergeChromosomes.internal ***
#
# DESCRIPTION:
#	subfunction of segragateChromosomes.internal, merging multiple chromosomes into one
# 
# PARAMETERS:
# 	cross - object of class cross
# 	chromosomes - chromosomes to be merged
# 	name - name of merged chromosome
# 
# OUTPUT:
#	object of class cross
#
############################################################################################################
mergeChromosomes.internal <- function(cross, chromosomes, name, verbose=FALSE){
	if(verbose)cat("Merging chromosomes",chromosomes,"to form chromosome",name,"names:",names(cross$geno),"\n")
	geno <- cross$geno
	markerNames <- NULL
	for(j in chromosomes){
		if(j!=name) markerNames <- c(markerNames, colnames(geno[[j]]$data))
	}
	for(k in markerNames) cross <- movemarker(cross, k, name)
	invisible(cross)
}

############################################################################################################
#									*** switchChromosomes.internal ***
#
# DESCRIPTION:
#	switching two chromosomes of cross object
# 
# cross - object of R/qtl cross type
# chr1, chr2 - numbers of chromosomes to be switched (1,2) == (2,1)
#
############################################################################################################
switchChromosomes.internal <- function(cross, chr1, chr2){
	cat(chr1,chr2,"\n")
	if(chr1!=chr2){
		geno <- cross$geno
		cross$geno[[chr1]] <- geno[[chr2]] 
		cross$geno[[chr2]] <- geno[[chr1]]
		cross <- est.rf(cross)
	}
	invisible(cross)
}


############################################################################################################
#									*** reduceChromosomesNumber ***
#
# DESCRIPTION:
#	Function to remove chromosomes from cross object. Those can specified in three ways described below.
# 
# PARAMETERS:
# 	cross - object of class cross
# 	numberOfChromosomes - how many chromosomes should stay (remove all but 1:numberOfChromosomes)
#	verbose - be verbose
# 
# OUTPUT:
#	object of class cross
#
############################################################################################################
reduceChromosomesNumber <- function(cross, numberOfChromosomes,verbose=FALSE){
	if(is.null(cross)&&!(any(class(cross)=="cross"))) stop("Not a cross object!\n")
	if(!(missing(numberOfChromosomes))){
		if(numberOfChromosomes<length(cross$geno)){
			for(i in length(cross$geno):(numberOfChromosomes+1)){
				cross <- removeChromosomesSub.internal(cross,i,verbose)
			}
		}
	}else{
		stop("You have to provide one of following: numberOfChromosomes, chromosomes or minLength")
	}
	invisible(cross)
}

############################################################################################################
#									*** removeChromosomes ***
#
# DESCRIPTION:
#	Function to remove chromosomes from cross object. Those can specified in three ways described below.
# 
# PARAMETERS:
# 	numberOfChromosomes - how many chromosomes should stay (remove all but 1:numberOfChromosomes)
# 	chromosomesToBeRmv - explicitly provide functions with NAMES of chromosomes to be removed
#	verbose - be verbose
# 
# OUTPUT:
#	object of class cross
#
############################################################################################################
removeChromosomes <- function(cross, chromosomesToBeRmv, verbose=FALSE){
	if(is.null(cross)&&!(any(class(cross)=="cross"))) stop("Not a cross object!\n")
	if(!(missing(chromosomesToBeRmv))){
		for(i in chromosomesToBeRmv){
			if(!(i%in%names(cross$geno))){
				stop("There is no chromosome called ",i,"\n")
			}else{
				cross <- removeChromosomesSub.internal(cross,i,verbose)
			}
		}
	}else{
		stop("You have to provide one of following: numberOfChromosomes, chromosomes or minLength")
	}
	invisible(cross)
}

############################################################################################################
#									*** removeTooSmallChromosomes ***
#
# DESCRIPTION:
#	Function to remove chromosomes from cross object. Those can specified in three ways described below.
# 
# PARAMETERS:
# 	cross - object of class cross
#	verbose - be verbose
# 	minNrOfMarkers - specify minimal number of markers chromosome is allowed to have (remove all that have
#					 less markers than that)
# 
# OUTPUT:
#	object of class cross
#
############################################################################################################
removeTooSmallChromosomes <- function(cross, minNrOfMarkers, verbose=FALSE){
	if(is.null(cross)&&!(any(class(cross)=="cross"))) stop("Not a cross object!\n")
	if(!(missing(minNrOfMarkers))){
		if(length(cross$geno)>1){
			if(length(cross$geno[[1]]$map)<minNrOfMarkers) minNrOfMarkers <- length(cross$geno[[1]]$map)-1
			for(i in length(cross$geno):1){
				if(length(cross$geno[[i]]$map)<minNrOfMarkers){
					cross <- removeChromosomesSub.internal(cross,i,verbose)
				}
			}
		}
	}else{
		stop("You have to provide one of following: numberOfChromosomes, chromosomes or minLength")
	}
	invisible(cross)
}

############################################################################################################
#									*** removeChromosomesSub.internal ***
#
# DESCRIPTION:
#	subfunction of removeChromosomes.internal, removing from given cross object specified chromosome
# 
# PARAMETERS:
# 	cross - object of class cross
# 	chr - chromosome to be removed (number or name)
# 
# OUTPUT:
#	object of class cross
#
############################################################################################################
removeChromosomesSub.internal <- function(cross, chr,verbose=FALSE){
	if(verbose)cat("removing chromosome:",chr," markers:",names(cross$geno[[chr]]$map),"\n")
	cross$rmv <- cbind(cross$rmv,cross$geno[[chr]]$data)
	cross <- drop.markers(cross, names(cross$geno[[chr]]$map))
	invisible(cross)
}

############################################################################################################
#									*** smooth.geno ***
#
# DESCRIPTION:
#	checking if fitted normal distributions do not overlap
# 
# PARAMETERS:
# 	offspring - currently processed row
# 	EM - output of normalmixEM function
# 	overlapInd - how many individuals are allowed to be overlapping between distributions
# 
# OUTPUT:
#	boolean
#
############################################################################################################
smooth.geno <- function(cross,windowSize=1,chr,population,map=c("genetic","physical"),verbose=FALSE){
	if(!any(class(cross) == "cross")) stop("Input should have class \"cross\".")
  if(!(missing(population))){ 
    map <- checkParameters.internal(map,c("genetic","physical"),"map")
    if(map=="genetic"){
      matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$genetic))
      if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
      if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
        population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
        if(verbose) cat(nrow(population$offspring$genotypes$real)-length(matchingMarkers),"markers were removed due to name mismatch\n")
      }
      oldMarkers <- pull.geno(cross)[,rownames(population$offspring$genotypes$real)]
    }else{
      matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$physical))
      if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
      if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
        population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
        if(verbose) cat(nrow(population$offspring$genotypes$real)-length(matchingMarkers),"markers were removed due to name mismatch\n")
      }
      oldMarkers <- pull.geno(cross)[,rownames(population$offspring$genotypes$real)]
    }
  }else{
    oldMarkers <- NULL
  }
	n.ind <- nind(cross)
  #cross <- fill.geno(cross)
  if(missing(chr)) chr <- 1:nchr(cross)
  cross_geno  <- vector(length(chr),mode="list")
  for(i in 1:length(chr)){
  if(verbose)cat("--- chr",i,"----\n")
    cross_geno[[i]]  <- smooth.genoSub.internal(cross$geno[[chr[i]]],windowSize,oldMarkers,verbose)
  }
	for(i in 1:length(chr)){
		cross$geno[[chr[i]]]$data <- cross_geno[[i]]$data
	}
  cross <- fill.geno(cross)
  if(verbose)cat("running est.rf\n")
	cross <- est.rf(cross)
  if(verbose)cat("running est.map\n")
	cross <- recalculateMap.internal(cross)
	invisible(cross)
}

############################################################################################################
#									*** smooth.genoSub.internal ***
#
# DESCRIPTION:
#	checking if fitted normal distributions do not overlap
# 
# PARAMETERS:
# 	offspring - currently processed row
# 	EM - output of normalmixEM function
# 	overlapInd - how many individuals are allowed to be overlapping between distributions
# 
# OUTPUT:
#	boolean
#
############################################################################################################
smooth.genoSub.internal <- function(geno,windowSize,oldMarkers,verbose){
  if(!is.null(dim(geno$data))){
    if(ncol(geno$data)>windowSize){
    old_genotype <- geno$data
    old_genotype[which(is.na(old_genotype))]<- 0
    genotype <- old_genotype
    genotype <- t(apply(genotype,1,smooth.genoRow.internal,windowSize))
    if(verbose) cat("changed",sum(genotype!=old_genotype)/length(genotype)*100,"% values because of genotyping error\n")
    if(any(colnames(genotype)%in%colnames(oldMarkers))){
      markersToBeUnchanged <- colnames(genotype)[which(colnames(genotype)%in%colnames(oldMarkers))]
      #print(markersToBeUnchanged)
      #print(genotype[1:10,1:10])
      genotype[,markersToBeUnchanged] <- oldMarkers[,markersToBeUnchanged]
      #print(genotype[1:10,1:10])
    }
    geno$data <- genotype
    }
  }
	invisible(geno)
}

############################################################################################################
#									*** smooth.genoRow.internal ***
#
# DESCRIPTION:
#	checking if fitted normal distributions do not overlap
# 
# PARAMETERS:
# 	offspring - currently processed row
# 	EM - output of normalmixEM function
# 	overlapInd - how many individuals are allowed to be overlapping between distributions
# 
# OUTPUT:
#	boolean
#
############################################################################################################
smooth.genoRow.internal <- function(genoRow,windowSize){
  if(length(table(genoRow))>1){
    wrongMarkers <- NULL
    if(length(genoRow)>windowSize+2){
      if(any(genoRow[1:windowSize]!=genoRow[windowSize+1])&&any(genoRow[1:windowSize]!=genoRow[windowSize+2])&&(genoRow[windowSize+1]==genoRow[windowSize+2])){
        wrongMarkers<-c(wrongMarkers,1:windowSize)
      }
      for(i in 2:(length(genoRow)-windowSize-1)){
        if(any(genoRow[i:(i+windowSize)]!=genoRow[i-1])&&any(genoRow[i:(i+windowSize)]!=genoRow[i+windowSize+1])&&(genoRow[i-1]==genoRow[i+windowSize+1])){
          wrongMarkers<-c(wrongMarkers,i:(i+windowSize))
        }
      }
      if(any(genoRow[(length(genoRow)-windowSize):(length(genoRow))]!=genoRow[(length(genoRow)-windowSize-1)])&&any(genoRow[(length(genoRow)-windowSize):(length(genoRow))]!=genoRow[(length(genoRow)-windowSize-2)])&&(genoRow[(length(genoRow)-windowSize-1)]==genoRow[(length(genoRow)-windowSize-2)])){
        wrongMarkers<-c(wrongMarkers,(length(genoRow)-windowSize):(length(genoRow)))
      }
    }
    if(length(wrongMarkers)>0){
      genoRow[wrongMarkers]<-NA
    }
  }
  invisible(genoRow)
}

recalculateMap.internal <- function(cross,...){
	newmap <- est.map(cross,offset=0,...)
	cross2 <- replace.map(cross, newmap)
	invisible(cross2)
}
