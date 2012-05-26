#
# smooth.geno.r
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Nov, 2011
# Contains: smooth.geno
#

# smooth.geno
#
# DESCRIPTION:
#  Checking if fitted normal distributions do not overlap
# PARAMETERS:
#  - offspring - currently processed row
#  - EM - output of normalmixEM function
#  - overlapInd - how many individuals are allowed to be overlapping between distributions
# OUTPUT:
#  Boolean
#
smooth.geno <- function(cross,windowSize=1,chr,population,map=c("genetic","physical"),fix.positions=FALSE,verbose=FALSE,...){
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
      oldPositions <- population$maps$genetic
      oldMarkers <- pull.geno(cross)[,rownames(population$offspring$genotypes$real)]
    }else{
      matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$physical))
      if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
      if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
        population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
        if(verbose) cat(nrow(population$offspring$genotypes$real)-length(matchingMarkers),"markers were removed due to name mismatch\n")
      }
      oldPositions <- population$maps$physical
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
    cross <- recalculateMap.internal(cross,...)
  if(!missing(population) && fix.positions){
        for(i in 1:nchr(cross)){
            toBeFixed <- names(cross$geno[[i]]$map)[which(names(cross$geno[[i]]$map)%in%rownames(oldPositions))]
            print(cross$geno[[i]]$map[toBeFixed])
            cross$geno[[i]]$map[toBeFixed] <- oldPositions[toBeFixed,2]
            print("+++++++++++++++")
            print(cross$geno[[i]]$map[toBeFixed])
        }
    }
  invisible(cross)
}

############################################################################################################
#                  *** smooth.genoSub.internal ***
#
# DESCRIPTION:
#  checking if fitted normal distributions do not overlap
# 
# PARAMETERS:
#   offspring - currently processed row
#   EM - output of normalmixEM function
#   overlapInd - how many individuals are allowed to be overlapping between distributions
# 
# OUTPUT:
#  boolean
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
#                  *** smooth.genoRow.internal ***
#
# DESCRIPTION:
#  checking if fitted normal distributions do not overlap
# 
# PARAMETERS:
#   offspring - currently processed row
#   EM - output of normalmixEM function
#   overlapInd - how many individuals are allowed to be overlapping between distributions
# 
# OUTPUT:
#  boolean
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
