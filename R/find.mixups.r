############################################################################################################
#
# find.mixups.R
#
# Copyright (c) 2012, Konrad Zych
#
# Modified by Danny Arends
# 
# first written January 2012
# last modified January 2012
# last modified in version: 1.0.0
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
# Contains: cross.saturate, rearrangeMarkers,
#             bestCorelated.internal, map2mapCorrelationMatrix.internal, map2mapImage
#
############################################################################################################

###########################################################################################################
#                                    *** find.mixups ***
#
# DESCRIPTION:
# 	saturate existing genetic map adding markers derived from gene expression
# OUTPUT:
#	object of class cross
############################################################################################################
find.mixups <- function(population,map=c("genetic","physical"),n.qtls=50,threshold=15,verbose=FALSE){
  if(missing(population)) stop("Please provide a population object\n")
  check.population(population)
  if(is.null(population$offspring$genotypes$real)){
    stop("No original genotypes in population$offspring$genotypes$real, load them in using intoPopulation function\n")
  }
  if((n.qtls*10)<nrow(population$offspring$phenotypes)){
    n.selectedPhenotypes <- n.qtls*10
  }else if(n.qtls<nrow(population$offspring$phenotypes)){
    n.selectedPhenotypes <- n.qtls
  }else{
    stop("n.qtls parameter is too high! Please choose at maximum number of the phenotypes/10.\n")
  }
  if(threshold>100){
    warning("Too high threshold, none of the markers will pass it.\n")
  }
  #selectedphenotypes <- round(runif(n.selectedPhenotypes,1,nrow(population$offspring$phenotypes)))
  map <- checkParameters.internal(map,c("genetic","physical"),"map")
  if(map=="genetic"){
    matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$genetic))
    if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
    if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
      population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
      population$maps$genetic <- population$maps$genetic[rownames(population$offspring$genotypes$real),]
      if(verbose) cat(nrow(population$offspring$genotypes$real)-length(matchingMarkers),"markers were removed due to name mismatch\n")
    }
    population10pheno <- population
    population10pheno$offspring$phenotypes <- population10pheno$offspring$phenotypes[1:10,]
    aa <- tempfile()
    sink(aa)
    returncross <- genotypesToCross.internal(population10pheno,"real","map_genetic")
    returncross$pheno <- t(population$offspring$phenotypes)
    sink()
    file.remove(aa)
  }else{
    matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$physical))
    if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
    if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
      population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
      population$maps$physical <- population$maps$physical[rownames(population$offspring$genotypes$real),]
      if(verbose) cat(nrow(population$offspring$genotypes$real)-length(matchingMarkers),"markers were removed due to name mismatch\n")
    }
    #for faster creation of cross
    population10pheno <- population
    population10pheno$offspring$phenotypes <- population10pheno$offspring$phenotypes[1:10,]
    aa <- tempfile()
    sink(aa)
    returncross <- genotypesToCross.internal(population10pheno,"real","map_physical")
    returncross$pheno <- t(population$offspring$phenotypes)
    sink()
    file.remove(aa)
  }
  returncross <- calc.genoprob(returncross)
  i <- round(runif(1,1,nrow(population$offspring$phenotypes)))
  qtls_found <- 0
  qtls <- NULL
  markers <- rownames(population$offspring$phenotypes)
  scores <- vector(mode="numeric",length=ncol(population$offspring$phenotypes))
  names(scores) <- colnames(population$offspring$phenotypes)
  while(qtls_found<n.qtls){
    cur_phenotype <- matrix(scanone(returncross,pheno.col=i,method="hk")[,3],1,nrow(population$offspring$genotypes$real))
    cur_peaks <- getpeaks.internal(abs(cur_phenotype),threshold)
    if(any(cur_peaks==2)){
      peakLocations <- which(cur_peaks==2)
      qtls_found <- qtls_found + length(peakLocations)
      if(verbose) cat(qtls_found,"qtls found\n")
      old_names <- names(qtls)
      qtls <- c(qtls,rep(i,length(peakLocations)))
      names(qtls) <- c(old_names,peakLocations)
    }
    i <- round(runif(1,1,nrow(population$offspring$phenotypes)))
  }
  for(j in 1:qtls_found){
    group_a <- which(population$offspring$genotypes$real[as.numeric(names(qtls))[j],]==1)
    group_b <- which(population$offspring$genotypes$real[as.numeric(names(qtls))[j],]==2)
    scores <- scoreMixups.internal(group_a,group_b,scores,qtls_found,population$offspring$phenotypes[qtls[j],])
  }
  if(verbose){
    if(any(scores>50)){
      flagged <- which(scores>50)
      cat("Found",length(flagged),"possible mix-ups:\n")
      for(flag in flagged){      
        cat(names(scores)[flag],":",scores[flag],"% flagged\n")
      }
    }
  }
  invisible(scores)
}

scoreMixups.internal <- function(group_a,group_b,scores,qtls_found,curRow){
  meanGroupA <- mean(curRow[group_a])
  meanGroupB <- mean(curRow[group_b])
  rowMean <- mean(meanGroupA,meanGroupB)
  increase <- 1#/qtls_found*100
  if(meanGroupA>meanGroupB){
    if(any(curRow[group_a]<rowMean)){
      positions <- names(curRow[group_a])[which(curRow[group_a]<rowMean)]
      scores[positions] <- scores[positions]+increase
    }
    if(any(curRow[group_b]>rowMean)){
      positions <- names(curRow[group_b])[which(curRow[group_b]>rowMean)]
      scores[positions] <- scores[positions]+increase
    }
  }else{
    if(any(curRow[group_a]>rowMean)){
      positions <- names(curRow[group_a])[which(curRow[group_a]>rowMean)]
      scores[positions] <- scores[positions]+increase
    }
    if(any(curRow[group_b]<rowMean)){
      positions <- names(curRow[group_b])[which(curRow[group_b]<rowMean)]
      scores[positions] <- scores[positions]+increase
    }
  }
  invisible(scores)
}
