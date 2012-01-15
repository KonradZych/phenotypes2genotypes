############################################################################################################
#
# cross.denovo.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified December 2011
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
#      switchChromosomes.internal, removeChromosomes.internal, removeChromosomesSub.internal,
#
############################################################################################################


############################################################################################################
#                                            *** cross.denovo ***
#
# DESCRIPTION:
#  ordering chromosomes using genetic/physical map and majority rule
# OUTPUT:
#  an object of class cross
############################################################################################################
cross.denovo <- function(population, cross, n.chr, map=c("none","genetic","physical"), comparisonMethod = c(sumMajorityCorrelation,majorityCorrelation,meanCorrelation,majorityOfMarkers), 
assignFunction=c(assignMaximumNoConflicts,assignMaximum), reOrder=TRUE, use.orderMarkers=FALSE, verbose=FALSE, debugMode=0){
  #checks
  if(missing(population)) stop("provide population object\n")
  check.population(population)
  if(missing(cross)) cross <- cross.denovo.internal(population,n.chr,verbose=TRUE,debugMode=2)
  if(length(cross$geno)<=1) stop("Selected cross object contains too little chromosomes to proceed.")
  map <- checkParameters.internal(map,c("none","genetic","physical"),"map")
  comparisonMethod <- checkParameters.internal(comparisonMethod,c(sumMajorityCorrelation,majorityCorrelation,meanCorrelation,majorityOfMarkers),"comparisonMethod")
  assignFunction <- checkParameters.internal(assignFunction,c(assignMaximumNoConflicts,assignMaximum),"assignFunction")
  if(verbose && debugMode==1) cat("cross.denovo starting withour errors in checkpoints.\n")

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
  if(map=="genetic"){
    originalMap <- population$maps$genetic
  }
  if(map=="physical"){
    originalMap <- population$maps$physical
  }
  chromToChromArray <- comparisonMethod(cross, originalMap, population)
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
  if(use.orderMarkers){
    if(verbose)cat("Ordering markers inside the cross object\n")
    s0 <- proc.time()
    nmarkersPerChr <- nmar(cross2)
    nChr <- length(nmarkersPerChr)
    for(i in 1:nChr){
      cross <- orderMarkers(cross2,use.ripple=F,chr=i,verb=T)
      e1 <- proc.time()
      if(i<nChr) te <- ((e1-s0)[3]/sum(nmarkersPerChr[1:i]))*sum(nmarkersPerChr[(i+1):nChr])
      else te <- 0
      if(verbose) cat("Done ordering chromosome",i,"/",nChr,"Time remaining:",te,"seconds.\n")
    }
    e0 <- proc.time()
    if(verbose && debugMode==2)cat("Ordering markers inside the cross object done in:",(e0-s0)[3],"seconds.\n")
  }
    invisible(cross)
  }
}


############################################################################################################
#                  *** assignMaximum ***
#
# DESCRIPTION:
#  function returning for rows or cols - which element is having max value
# OUTPUT:
#  vector of numerics
############################################################################################################
assignMaximum <- function(x, use=2){
  apply(x,use,which.max)
}

############################################################################################################
#                  *** assignMaximumNoConflicts ***
#
# DESCRIPTION:
#  function return
# OUTPUT:
#  vector of numerics
############################################################################################################
assignMaximumNoConflicts <- function(x, use=2){
  assignment <- assignMaximum(x,use)
  notYetAssigned <- as.numeric(names(assignment)[which(!(names(assignment)%in%assignment))])
  while(any(duplicated(assignment))){
    duplicated_ones <- assignment[(duplicated(assignment))]
    for(duplication in duplicated_ones){
      need_to_decide <- which(assignment == duplication)
      best_fitting <- names(which.max(apply(x,use,max)[need_to_decide]))
      need_to_decide <- need_to_decide[-which(names(need_to_decide)==best_fitting)]
      #### we need to work on the next line, probably some while loop here, because if you call function
      #### few times(I mean all the functions on the way not only this one, it's getting better and better!:P
      assignment[as.numeric(names(need_to_decide))] <- notYetAssigned[1]
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
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
#
############################################################################################################
majorityCorrelation <- function(cross, originalMap, population, verbose=FALSE){
  nrOfChromosomesInCross <- nchr(cross)
  genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross, population, FALSE)
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
#   Subfunction of assignChromosomes, for every chromosome in cross for every marker checks the marker it is
#   having highest correlation with. Checks on which chromosome this marker is placed in old map. For each of
#   new chromosomes one or more of chromosomes from old map will be represented. Function sums correlations for
#   each pair of those and for every new chromosomes assigns old chromosome with highest cumulative cor.
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
#
############################################################################################################
sumMajorityCorrelation <- function(cross, originalMap, population, verbose=FALSE){
  nrOfChromosomesInCross <- nchr(cross)
  genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross, population, FALSE)
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

############################################################################################################
#                                         *** correlationRule ***
#
# DESCRIPTION:
#  Subfunction of assignChromosomes, assigning chromosome from new map to old ones using mean corelation\
#  between their markers.
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
#
############################################################################################################
meanCorrelation <- function(cross,originalMap,population,verbose=FALSE){
  nrOfChromosomesInCross <- nchr(cross)
  genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross, population, FALSE)
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
majorityOfMarkers <- function(cross,originalMap,population,verbose=FALSE){
  nrOfChromosomesInCross <- nchr(cross)
  chromToChromArray <- matrix(0,length(unique(originalMap[,1])),nrOfChromosomesInCross)
  rownames(chromToChromArray) <- unique(originalMap[,1])
  colnames(chromToChromArray) <- 1:nrOfChromosomesInCross
  for(i in 1:nrOfChromosomesInCross){
    markersFromCurNewChrom <-colnames(cross$geno[[i]]$data)
    for(j in unique(originalMap[,1])){
        markersFromCurOldChrom <- rownames(originalMap)[which(originalMap[,1]==j)]
        chromToChromArray[j,i] <- sum(markersFromCurOldChrom%in%markersFromCurNewChrom)
    }
  }
  return(chromToChromArray)
}


############################################################################################################
#                  *** cross.denovo.internal ***
#
# DESCRIPTION:
#   function to create new map and save it in cross object
# 
# PARAMETERS:
#   population - object of class population
#   orde - object of class population
#   n.chr - expected number of linkage groups
#   use - expected number of linkage groups
#  verbose - be verbose
#
# OUTPUT:
#  an object of class cross
#
#
############################################################################################################
cross.denovo.internal<- function(population,  n.chr,  use=c("geno","rf"), verbose=FALSE, debugMode=0){
  if(missing(n.chr)) stop("n.chr in an obligatory parameter")
  if(missing(population)) stop("no population object provided")
  use <- checkParameters.internal(use,c("geno","rf"),"use")
  check.population(population)
  if(is.null(population$offspring$genotypes$simulated)){
    stop("no simulated genotypes in population object, first use findBiomarkers!\n")
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
