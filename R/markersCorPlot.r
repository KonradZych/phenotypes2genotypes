############################################################################################################
#
# markersCorPlot.R
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
# Contains: markersCorPlot, ascendingMaptoJigSawMap
#               getChrOffsets.internal, getMarkerOffsets, getMarkerOffsetsFromMap, 
#               getPopulationOffsets.internal, chromCorMatrix
#
############################################################################################################

############################################################################################################
#									*** markersCorPlot ***
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
markersCorPlot <- function(cross, population, map=c("genetic","physical"), cmBetween=25, comparisonMethod = c(sumMajorityCorrelation,majorityCorrelation,meanCorrelation), chr,verbose=TRUE){
  ### checks
  map <- defaultCheck.internal(map,"map",2,"genetic")
	if(map=="genetic"){
    originalMap <- population$maps$genetic
  }else{
    originalMap <- population$maps$physical
  }
  comparisonMethod <- defaultCheck.internal(comparisonMethod,"comparisonMethod",3,sumMajorityCorrelation)
  if(is.null(originalMap)) stop("no ",map," map provided!")  
  
  ### getting offsets for each chromosome on both maps
  offsets1 <- getPopulationOffsets.internal(population,originalMap,cmBetween)
  n.originalChrom <- length(offsets1)-1
  offsets2 <- getChrOffsets.internal(cross,cmBetween)
  n.newChrom <- length(offsets2)-1
  
  if(n.originalChrom<n.newChrom){
	offsets1 <- c(offsets1,rep(0,(n.newChrom- n.originalChrom)))
  }else if( n.originalChrom>n.newChrom){
	offsets2 <- c(offsets2,rep(0,( n.originalChrom-n.newChrom)))
  }
    
  ### global offsets
  global_offset <- NULL
  for(x in 1:length(offsets1)){
    global_offset <- c(global_offset,max(offsets1[x],offsets2[x]))
  }
  
  ### positions of markers (absolute - with offsets)
  mloc_original <- getMarkerOffsets(cross,global_offset[1:n.newChrom],cmBetween)
  mloc_o <- getMarkerOffsetsFromMap(originalMap,global_offset[1:n.originalChrom],cmBetween)

  ### limits of plot

  
  ### summary offsets
  sum_gl_off <- NULL
  for(x in 1:length(global_offset)){
    sum_gl_off <- c(sum_gl_off,sum(global_offset[1:x]))
  }
  if(missing(chr)){
    m_max <- max(mloc_o,mloc_original)
    m_min <- min(mloc_o,mloc_original)
  }else{
    m_max <- sum_gl_off[max(chr)+1]
    m_min <- sum_gl_off[min(chr)]
  }
  cat("-----1-----\n")
  ### preparing chrom to chrom cor matrix for use in the background
  genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross, population, FALSE)
  chromToChromMatrix <- comparisonMethod(cross,originalMap,population)
  maximum <- max(chromToChromMatrix)
  cat("-----2-----\n")
  
  ### setting plot canvas
  plot(c(m_min,m_max),c(m_min,m_max),type='n',xlab="Original map",ylab="New map",main="Comparison of genetic maps", xaxt="n", yaxt="n")
  ### background
  for(i in 1:(n.originalChrom)){
    for(j in 1:(n.newChrom)){
        cur_col <- (maximum-(chromToChromMatrix[i,j]))/maximum
        rect(sum_gl_off[i],sum_gl_off[j],sum_gl_off[i+1],sum_gl_off[j+1],lty=0,col=rgb(cur_col,cur_col,cur_col))
    }
  }
  ### markers on new map
  points(cbind(mloc_o,mloc_o),pch=21,col="red",cex=1.5,lwd=4)
  ### markers on original map
  points(cbind(mloc_original,mloc_original),pch=20,col="green",cex=1.5,lwd=2)
  ### gris
  abline(v=sum_gl_off[-length(sum_gl_off)],lty=2)
  abline(h=sum_gl_off[-length(sum_gl_off)],lty=2)
  ### chromosome labels and tics
  labelsPos <- vector(mode="numeric",length(sum_gl_off)-1)
	for(i in 1:length(sum_gl_off)-1){
		labelsPos[i] <- (sum_gl_off[i] + sum_gl_off[i+1])/2
	}
	axis(1, at = sum_gl_off[-length(sum_gl_off)],labels = FALSE)
	axis(1, at = labelsPos[1:n.originalChrom],labels = names(table(originalMap[,1])), lwd = 0, tick = FALSE)
	axis(2, at = sum_gl_off[-length(sum_gl_off)],labels = FALSE)
	axis(2, at = labelsPos[1:n.newChrom],labels = chrnames(cross), lwd = 0, tick = FALSE)
  invisible(chromToChromMatrix)
}

############################################################################################################
#									*** markersCorPlot ***
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
getChrOffsets.internal <- function(cross, cmBetween){
  offsets <- unlist(lapply(pull.map(cross),max))
  offsets <- offsets+cmBetween
  offsets <-c(0,offsets)
  offsets
}


#From IQTL by Danny Arends, SHOULD NOT MODIFY
getMarkerOffsets <- function(cross, offsets, cmBetween=25){
  if(missing(offsets))offsets <- getChrOffsets.internal(cross,cmBetween)
  cnt <- 1
  myoffsets <- NULL
  for(x in nmar(cross)){
    myoffsets <- c(myoffsets,rep(sum(offsets[1:cnt]),x))
    cnt <- cnt+1
  }

  mlocations <- myoffsets + as.numeric(unlist(pull.map(cross)))
  mlocations
}

############################################################################################################
#									*** markersCorPlot ***
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
getMarkerOffsetsFromMap <- function(map, offsets, cmBetween=25){
  cnt <- 1
  myoffsets <- NULL
  for(x in table(map[,1])){
    myoffsets <- c(myoffsets,rep(sum(offsets[1:cnt]),x))
    cnt <- cnt+1
  }

  mlocations <- myoffsets + as.numeric(map[,2])
  mlocations
}


############################################################################################################
#									*** markersCorPlot ***
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
ascendingMaptoJigSawMap <- function(mapToProcess,verbose=FALSE){
  if(is.null(mapToProcess)) stop("No map provided!")
  for(x in unique(mapToProcess[,1])){
    if(verbose) cat("Processing chromosome:",x,"\n")
    offsetOfCurrentChromosome <- min(mapToProcess[which(mapToProcess[,1]==x),2])
    mapToProcess[which(mapToProcess[,1]==x),2] <- mapToProcess[which(mapToProcess[,1]==x),2]-offsetOfCurrentChromosome
  }
  invisible(mapToProcess)
}

############################################################################################################
#									*** markersCorPlot ***
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
getPopulationOffsets.internal <- function(population, originalMap, cmBetween){
  minima <- NULL
  for(x in unique(originalMap[,1])){
    minima <- c(minima,max(originalMap[which(originalMap[,1]==x),2]))
  }
  minima <- minima + cmBetween
  minima <- c(0,minima)
  invisible(minima)
}

############################################################################################################
#									*** markersCorPlot ***
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
chromCorMatrix <- function(cross,population,map=c("genetic","physical"),show=c(max,mean),verbose=FALSE){
  map <- defaultCheck.internal(map,"map",2,"genetic")
	if(map=="genetic"){
    old_map <- population$maps$genetic
  }else{
    old_map <- population$maps$physical
  }
  if(is.null(old_map)) stop("no ",map," map provided!")  
  result <- matrix(0,nchr(cross),length(unique(old_map[,1])))
  genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross,population,verbose)
  s <- proc.time()
  for(i in 1:nchr(cross)){
      markersfromnewmap <- colnames(cross$geno[[i]]$data)
      for(j in unique(old_map[,1])){
         rownamesOfSomthing <- rownames(old_map)[which(old_map[,1]==j)]
         markersfromoldmap <- rownames(population$offspring$genotypes$real[rownamesOfSomthing,])
         result[i,j]<- show(abs(genotypesCorelationMatrix[markersfromoldmap,markersfromnewmap]))
      }
      e <- proc.time()
      cat("Counting correlation matrix:",round(i/nchr(cross)*100),"% done estimated time remaining:",((e-s)[3]/i)*(nchr(cross)-i),"s\n")
  }
  invisible(result)
}