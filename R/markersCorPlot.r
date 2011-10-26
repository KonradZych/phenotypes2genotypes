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
markersCorPlot <- function(cross, population, map=c("genetic","physical"), cmBetween=25, show=c(max,mean),verbose=TRUE){
  map <- defaultCheck.internal(map,"map",2,"genetic")
	if(map=="genetic"){
    cur_map <- population$maps$genetic
  }else{
    cur_map <- population$maps$physical
  }
  if(is.null(cur_map)) stop("no ",map," map provided!")  
  offsets1 <- getPopulationOffsets.internal(population,cur_map,cmBetween)
  offsets2 <- getChrOffsets.internal(cross,cmBetween)

  global_offset <- NULL
  for(x in 1:length(offsets1)){
    global_offset <- c(global_offset,max(offsets1[x],offsets2[x]))
  }
  
  sum_gl_off <- NULL
  for(x in 1:length(global_offset)){
    sum_gl_off <- c(sum_gl_off,sum(global_offset[1:x]))
  }

  mloc_original <- getMarkerOffsets(cross,global_offset,cmBetween)
  mloc_o <- getMarkerOffsetsFromMap(cur_map,global_offset,cmBetween)

  m_max <- max(mloc_o,mloc_original)
  m_min <- min(mloc_o,mloc_original)

  back <- chromCorMatrix(cross,population,map,show,verbose)
  plot(c(m_min,m_max),c(m_min,m_max),type='n',xlab="Old map (cM)",ylab="New map (cM)",main="Plot comparison")
  for(i in 1:(length(sum_gl_off)-1)){
    for(j in 1:(length(sum_gl_off)-1)){
        cur_col <- 1-(back[i,j])
        rect(sum_gl_off[i],sum_gl_off[j],sum_gl_off[i+1],sum_gl_off[j+1],lty=0,col=rgb(cur_col,cur_col,cur_col))
    }
  }
  points(cbind(mloc_o,mloc_o),pch=21,col="red",lwd=4)
  points(cbind(mloc_original,mloc_original),pch=20,col="green",lwd=2)
  abline(v=sum_gl_off[-length(sum_gl_off)],lty=2)
  abline(h=sum_gl_off[-length(sum_gl_off)],lty=2)
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
getPopulationOffsets.internal <- function(population, cur_map, cmBetween){
  minima <- NULL
  for(x in unique(cur_map[,1])){
    minima <- c(minima,max(cur_map[which(cur_map[,1]==x),2]))
  }
  minima <- minima + cmBetween
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