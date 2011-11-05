############################################################################################################
#
# postprocessing.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified September 2011
# last modified in version: 0.9.0
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
#									*** assignChromosomes ***
#
# DESCRIPTION:
#	ordering chromosomes using genetic/physical map and majority rule
# 
# PARAMETERS:
# 	cross - object of class cross, containing physical or genetic map
# 	map - which map should be used for comparison:
#			- genetic - genetic map from cross$maps$genetic
#			- physical - physical map from cross$maps$physical
#	verbose - be verbose
#
# OUTPUT:
#	object of class cross
#
############################################################################################################
assignChromosomes <- function(population, cross, n.chr, map=c("none","genetic","physical"), comparisonMethod = c(sumMajorityCorrelation,majorityCorrelation,meanCorrelation), assignFunction=c(assignMaximumNoConflicts,assignMaximum),reOrder=FALSE, verbose=FALSE, orderMarkersOnChrom=FALSE, debugMode=0){
  if(missing(cross)){
    cat("Cross object not found, will be created form population object\n")
    cross <- createNewMap(population,n.chr,verbose=TRUE,debugMode=2)
  }
  if(length(cross$geno)<=1) stop("selected cross object contains too little chromosomes to proceed")
  map <- defaultCheck.internal(map,"map",3,"none")
  comparisonMethod <- defaultCheck.internal(comparisonMethod,"comparisonMethod",3,sumMajorityCorrelation)
  assignFunction <- defaultCheck.internal(assignFunction,"assignFunction",2,assignMaximumNoConflicts)
  if(map=="none"){
    if(reOrder){
      return(cross)
    }else{
      assignment <- names(cross$geno)
      names(assignment) <- names(cross$geno)
      return(assignment)
    }
  }
  
  s1 <- proc.time()
  genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross, population, FALSE)

  if(map=="genetic"){
    originalMap <- population$maps$genetic
  }
  if(map=="physical"){
    originalMap <- population$maps$physical
  }
  chromToChromArray <- comparisonMethod(cross, originalMap, genotypesCorelationMatrix)
  e1 <- proc.time()
  if(verbose)cat("Calculating correlation matrix done in:",(e1-s1)[3],"seconds.\n")
  
  assignment <- assignFunction(chromToChromArray)
  
  if(!reOrder){
    if(verbose)cat("Returning new ordering vector.\n")
    invisible(assignment)
  }else{
    ordering <- assignedChrToMarkers(assignment,cross)
    if(verbose)cat("Applying new ordering to the cross object.\n")
    cross2 <- reorganizeMarkersWithin(cross,ordering)
    if(verbose)cat("Ordering markers inside the cross object\n")
    s1 <- proc.time()
    aa <- tempfile()
    sink(aa)
    if(orderMarkersOnChrom) cross <- orderMarkers(cross2,use.ripple=F,verb=T)
    sink()
    file.remove(aa)
    e1 <- proc.time()
    if(verbose && debugMode==2)cat("Saving data into cross object done in:",(e1-s1)[3],"seconds.\n")
    invisible(cross2)
  }
}



assignMaximum <- function(x, use=1){
  apply(x,use,which.max)
}

assignMaximumNoConflicts <- function(x, use=1){
  assignment <- assignMaximum(x,use)
  notYetAssigned <- as.numeric(names(assignment)[which(!(names(assignment)%in%assignment))])
  while(any(duplicated(assignment))){
    duplicated_ones <- assignment[(duplicated(assignment))]
    for(duplication in duplicated_ones){
      need_to_decide <- which(assignment == duplication)
      best_fitting <- names(which.max(apply(x,1,max)[need_to_decide]))
      need_to_decide <- need_to_decide[-which(names(need_to_decide)==best_fitting)]
      #### we need to work on the next line, probably some while loop here, because if you call function
      #### few times(I mean all the functyions on the way not aonly this one, it's getting better and better!:P
      assignment[as.numeric(names(need_to_decide))] <- notYetAssigned[1:length(need_to_decide)]
      notYetAssigned <- as.numeric(names(assignment)[which(!(names(assignment)%in%assignment))])
      }
    }
  invisible(assignment)
}

############################################################################################################
#                                           *** majorityRule ***
#
# DESCRIPTION:
#   Subfunction of assignChromosomes, for every chromosome in cross for every marker checks the marker it is
#   having highest correlation with. Checks on which chromosome this marker is placed in old map. For each of
#   new chromosomee, old chromosome with most markers with high correlation is assigned. 
# 
# PARAMETERS:
#  cross - object of class cross
#  originalMap - map from population object
#  genotypesCorelationMatrix - gene correlation matrix (from map2mapCorrelationMatrix function)
#  verbose - be verbose
# 
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
#
############################################################################################################
majorityCorrelation <- function(cross, originalMap, genotypesCorelationMatrix, verbose=FALSE){
  nrOfChromosomesInCross <- nchr(cross)
  chrMaxCorMarkers <- t(apply(abs(genotypesCorelationMatrix),2,function(r){c(originalMap[rownames(genotypesCorelationMatrix)[which.max(r)],1],max(r))}))
  ordering <- NULL
  chromToChromArray <- matrix(0,length(unique(originalMap[,1])),nrOfChromosomesInCross)
  rownames(chromToChromArray) <- unique(originalMap[,1])
  colnames(chromToChromArray) <- 1:nrOfChromosomesInCross
  for(i in 1:nrOfChromosomesInCross){
    markersFromCurChrom <-colnames(cross$geno[[i]]$data)
    #### USING ONLY MARKERS WITH COR HIGHER THAN 0.5
    markersCorWithCurChrom <- chrMaxCorMarkers[markersFromCurChrom,]
    if(length(markersFromCurChrom)>1){
        #markersHighlyCorWithCurChrom <- markersCorWithCurChrom[which(abs(markersCorWithCurChrom[,2])>0.5),]
        chromosomesCorWithCurrent <- table(markersCorWithCurChrom[,1])
        #bestCorChrom <- as.numeric(names(chromosomesCorWithCurrent)[which.max(chromosomesCorWithCurrent)])
        for(j in names(chromosomesCorWithCurrent)){
          chromToChromArray[j,i] <- chromosomesCorWithCurrent[j]
        }
    }else if(length(markersFromCurChrom)==1){
        #markersHighlyCorWithCurChrom <- markersCorWithCurChrom[which(abs(markersCorWithCurChrom[,2])>0.5),]
        chromToChromArray[markersCorWithCurChrom[1],i] <- markersCorWithCurChrom[2]
    }else{
         bestCorChrom <- NA
    }
    
  }
  invisible(chromToChromArray)
}



############################################################################################################
#                                          ** sumMajorityRule***
#
# DESCRIPTION:
# 	Subfunction of assignChromosomes, for every chromosome in cross for every marker checks the marker it is
#   having highest correlation with. Checks on which chromosome this marker is placed in old map. For each of
#   new chromosomes one or more of chromosomes from old map will be represented. Function sums correlations for
#   each pair of those and for every new chromosomes assigns old chromosome with highest cumulative cor.
# 
# PARAMETERS:
#  cross - object of class cross
#  originalMap - map from population object
#  genotypesCorelationMatrix - gene correlation matrix (from map2mapCorrelationMatrix function)
#  verbose - be verbose
# 
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
#
############################################################################################################
sumMajorityCorrelation <- function(cross, originalMap, genotypesCorelationMatrix, verbose=FALSE){
  nrOfChromosomesInCross <- nchr(cross)
  chrMaxCorMarkers <- t(apply(abs(genotypesCorelationMatrix),2,function(r){c(originalMap[rownames(genotypesCorelationMatrix)[which.max(r)],1],max(r))}))
  chromToChromArray <- matrix(0,length(unique(originalMap[,1])),nrOfChromosomesInCross)
  rownames(chromToChromArray) <- unique(originalMap[,1])
  colnames(chromToChromArray) <- 1:nrOfChromosomesInCross
  for(i in 1:nrOfChromosomesInCross){
    markersFromCurChrom <-colnames(cross$geno[[i]]$data)
    markersHighlyCorWithCurChrom <- chrMaxCorMarkers[markersFromCurChrom,]
    if(length(markersFromCurChrom)>1){
        correlatedChrom <- as.numeric(names(table(markersHighlyCorWithCurChrom[,1])))
        for(j in correlatedChrom){
              currentSum <- sum(markersHighlyCorWithCurChrom[which(markersHighlyCorWithCurChrom[,1]==j),2])
              chromToChromArray[j,i] <- currentSum
        }
    }else{
       chromToChromArray[markersHighlyCorWithCurChrom[1],i] <- markersHighlyCorWithCurChrom[2]
    }
  }
  invisible(chromToChromArray)
}

assignedChrToMarkers <- function(assignment,cross){
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
#                                         *** correlationRule ***
#
# DESCRIPTION:
#  Subfunction of assignChromosomes, assigning chromosome from new map to old ones using mean corelation\
#  between their markers.
# 
# PARAMETERS:
#  cross - object of class cross
#  originalMap - map from population object
#  genotypesCorelationMatrix - gene correlation matrix (from map2mapCorrelationMatrix function)
#  verbose - be verbose
# 
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
#
############################################################################################################
meanCorrelation <- function(cross,originalMap,genotypesCorelationMatrix,verbose=FALSE){
  nrOfChromosomesInCross <- nchr(cross)
  chromToChromArray <- matrix(0,length(unique(originalMap[,1])),nrOfChromosomesInCross)
  rownames(chromToChromArray) <- unique(originalMap[,1])
  colnames(chromToChromArray) <- 1:nrOfChromosomesInCross
  for(i in 1:nrOfChromosomesInCross){
    markersFromCurNewChrom <-colnames(cross$geno[[i]]$data)
    for(j in unique(originalMap[,1])){
        markersFromCurOldChrom <- rownames(originalMap)[which(originalMap[,1]==j)]
        currentCorMatrix <- abs(genotypesCorelationMatrix[markersFromCurOldChrom,markersFromCurNewChrom])
        chromToChromArray[j,i] <- sum(apply(currentCorMatrix,1,mean))
    }
  }
  return(chromToChromArray)
}

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
#									*** smoothGeno ***
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
smoothGeno <- function(cross,windowSize=1,chr,verbose=FALSE){
	if(!any(class(cross) == "cross")) stop("Input should have class \"cross\".")
	n.ind <- nind(cross)
  cross <- fill.geno(cross)
  if(missing(chr)) chr <- 1:nchr(cross)
  cross_geno  <- vector(length(chr),mode="list")
  for(i in 1:length(chr)){
  if(verbose)cat("--- chr",i,"----\n")
    cross_geno[[i]]  <- smoothGenoSub.internal(cross$geno[[chr[i]]],windowSize,verbose)
  }
	for(i in 1:length(chr)){
		cross$geno[[chr[i]]]$data <- cross_geno[[i]]$data
	}
  if(verbose)cat("running est.rf\n")
	cross <- est.rf(cross)
  if(verbose)cat("running est.map\n")
	cross <- recalculateMap.internal(cross)
	invisible(cross)
}

############################################################################################################
#									*** smoothGenoSub.internal ***
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
smoothGenoSub.internal <- function(geno,windowSize,verbose){
  if(!is.null(dim(geno$data))){
    if(ncol(geno$data)>windowSize){
    old_genotype <- geno$data
    genotype <- old_genotype
    genotype <- t(apply(genotype,1,smoothGenoRow.internal,windowSize))
    if(verbose) cat("changed",sum(genotype!=old_genotype)/length(genotype)*100,"% values because of genotyping error\n")
    geno$data <- genotype
    }
  }
	invisible(geno)
}

############################################################################################################
#									*** smoothGenoRow.internal ***
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
smoothGenoRow.internal <- function(genoRow,windowSize){
  if(length(genoRow)>windowSize+2){
    if(any(genoRow[1:windowSize]!=genoRow[windowSize+1])&&any(genoRow[1:windowSize]!=genoRow[windowSize+2])&&(genoRow[windowSize+1]==genoRow[windowSize+2])){
      genoRow[1:windowSize] <- genoRow[windowSize+1]
    }
    for(i in 2:(length(genoRow)-windowSize-1)){
      if(any(genoRow[i:(i+windowSize)]!=genoRow[i-1])&&any(genoRow[i:(i+windowSize)]!=genoRow[i+windowSize+1])&&(genoRow[i-1]==genoRow[i+windowSize+1])){
        genoRow[i:(i+windowSize)] <- genoRow[i-1]
      }
    }
    if(any(genoRow[(length(genoRow)-windowSize):(length(genoRow))]!=genoRow[(length(genoRow)-windowSize-1)])&&any(genoRow[(length(genoRow)-windowSize):(length(genoRow))]!=genoRow[(length(genoRow)-windowSize-2)])&&(genoRow[(length(genoRow)-windowSize-1)]==genoRow[(length(genoRow)-windowSize-2)])){
      genoRow[(length(genoRow)-windowSize):(length(genoRow))] <- genoRow[(length(genoRow)-windowSize-1)]
    }
  }
  invisible(genoRow)
}

recalculateMap.internal <- function(cross,...){
	newmap <- est.map(cross,offset=0,...)
	cross2 <- replace.map(cross, newmap)
	invisible(cross2)
}
