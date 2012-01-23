############################################################################################################
#
# cross.saturate.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
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
#                                    *** cross.saturate ***
#
# DESCRIPTION:
# 	saturate existing genetic map adding markers derived from gene expression
# OUTPUT:
#	object of class cross
############################################################################################################
cross.saturate <- function(population, cross, map=c("genetic","physical"), placeUsing=c("qtl","correlation"), threshold=3, use.orderMarkers=FALSE, verbose=FALSE, debugMode=0){
  if(missing(population)) stop("Please provide a population object\n")
  if(is.null(population$offspring$genotypes$real)){
    stop("No original genotypes in population$offspring$genotypes$real, load them in using intoPopulation\n")
  }
  check.population(population)
  if(!is.numeric(threshold)||is.na(threshold)) stop("Please provide correct threshold")
  if(threshold<=0){
    cat("WARNING: threshold too low, all possible markers will be selected\n")
  }else if(threshold>=5){
    cat("WARNING: threshold too high, few new markers will be selected\n")
  }
  map <- checkParameters.internal(map,c("genetic","physical"),"map")
  placeUsing <- checkParameters.internal(placeUsing,c("qtl","correlation"),"placeUsing")
 if(map=="genetic"){
    matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$genetic))
    if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
    if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
      population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
      population$maps$genetic <- population$maps$genetic[rownames(population$offspring$genotypes$real),]
      if(verbose) cat(nrow(population$offspring$genotypes$real)-length(matchingMarkers),"markers were removed due to name mismatch\n")
    }
    cur_map <- population$maps$genetic
  }else{
    matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$physical))
    if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
    if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
      population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
      population$maps$physical <- population$maps$physical[rownames(population$offspring$genotypes$real),]
      if(verbose) cat(nrow(population$offspring$genotypes$real)-length(matchingMarkers),"markers were removed due to name mismatch\n")
    }
    cur_map <- population$maps$physical
  }
  if(missing(cross)){
    if(is.null(population$offspring$genotypes$simulated)){
      stop("No genotype data in population$offspring$genotypes$simulated, run findBiomarkers first\n")
    }else{
      cat("No cross object provided, creating one using population object\n")
      	#*******SAVING CROSS OBJECT*******
      s1 <- proc.time()
      aa <- tempfile()
      sink(aa)
      cross <- genotypesToCross.internal(population,"simulated",verbose=verbose,debugMode=debugMode)
      sink()
      file.remove(aa)
      e1 <- proc.time()
      if(verbose && debugMode==2)cat("Saving data into cross object done in:",(e1-s1)[3],"seconds.\n")
    }
  }
 
  #*******ENRICHING ORIGINAL MAP*******
	s1 <- proc.time()
	cross <- rearrangeMarkers(cross,population,cur_map,threshold,placeUsing,addMarkers=TRUE,verbose=verbose)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Enrichment of original map done in:",(e1-s1)[3],"seconds.\n")
  
  #*******ORDERING NEW MAP*******
  if(use.orderMarkers){
    if(verbose)cat("Ordering markers inside the cross object\n")
    s1 <- proc.time()
    aa <- tempfile()
    sink(aa)
    cross <- orderMarkers(cross,use.ripple=F,verb=T)
    sink()
    file.remove(aa)
    e1 <- proc.time()
    if(verbose && debugMode==2)cat("Saving data into cross object done in:",(e1-s1)[3],"seconds.\n")
   }
	invisible(cross)
}


###########################################################################################################
#                                    *** rearrangeMarkers ***
#
# DESCRIPTION:
# 	ordering chromosomes using genetic/physical map and corelation rule
# OUTPUT:
#	object of class cross
############################################################################################################
rearrangeMarkers <- function(cross, population, cur_map, threshold=3, placeUsing,addMarkers=FALSE, verbose=FALSE){
  if(verbose) cat("old map contains",max(cur_map[,1]),"chromosomes\n")
  if(placeUsing=="qtl"){
    markersNewPostions <- bestQTL.internal(cross,population,threshold,verbose)
  }else{
    markersNewPostions <- bestCorelated.internal(cross,population,cur_map,threshold,verbose)
  }
  if(nrow(markersNewPostions) == 0){
    cat("selected",nrow(markersNewPostions),"markers with current corThreshold, there will be only markers from old map in the cross object\n")
  }else if(verbose){
    cat("selected",nrow(markersNewPostions),"markers for further analysis\n")
  }
	returncross <- cross
	returncross$geno <- vector(max(cur_map[,1]), mode="list")
	returncross$pheno <- pull.pheno(cross)
	if(verbose) cat("Reordering markers \n")  
	for(x in unique(cur_map[,1])){
		if(verbose) cat("- chr ",x," -\n")    
		oldnames <- rownames(cur_map)[which(cur_map[,1]==x)]
    oldpositions <- cur_map[oldnames,2]
    newnames <- rownames(markersNewPostions)[which(markersNewPostions[,1]==x)]
    newpositions <- markersNewPostions[newnames,2]
    toRmv <- NULL
    if(length(newnames)>0){
      for(i in 1:length(newpositions)){
        if(newpositions[i]%in%oldpositions){
          toRmv <- c(toRmv,i)
        }
      }
      if(length(toRmv)>0){
        newnames <- newnames[-toRmv]
        newpositions <- newpositions[-toRmv]
        }
    }
    if(verbose) cat("Selected:",length(newnames),"new and",length(oldnames),"original markers,",length(toRmv),"markers were removed\n") 
		if(addMarkers){
			returncross$geno[[x]]$data <- cbind(pull.geno(cross)[,newnames],t(population$offspring$genotypes$real[oldnames,]))
			newmap <- c(as.numeric(newpositions),oldpositions)
			names(newmap) <- c(newnames,oldnames)
      newmap <- sort(newmap)
      colnames(returncross$geno[[x]]$data) <- c(newnames,oldnames)
      returncross$geno[[x]]$data <- returncross$geno[[x]]$data[,names(newmap)]
		}else{
			returncross$geno[[x]]$data <- pull.geno(cross)[,newnames]
			newmap <- as.numeric(newpositions)
			names(newmap) <- newnames
      newmap <- sort(newmap)
      colnames(returncross$geno[[x]]$data) <- newnames
      returncross$geno[[x]]$data <- returncross$geno[[x]]$data[,names(newmap)]
		}
		returncross$geno[[x]]$map <- c(newmap)
	}
	names(returncross$geno) <- 1:length(returncross$geno)
	for(i in 1:length(returncross$geno)){
		class(returncross$geno[[i]]) <- "A"
	}
	invisible(returncross)
}


###########################################################################################################
#                                    *** bestCorelated.internal ***
#
# DESCRIPTION:
# 	subfunction of segragateChromosomes.internal, returns matrix showing for every reco map chromosome from 
#	which physicall map chromosome majority of markers comes
# OUTPUT:
#	vector with new ordering of chromosomes inside cross object
############################################################################################################
bestCorelated.internal <- function(cross,population, cur_map,corSDTreshold,verbose=FALSE){
  cormatrix <- map2mapCorrelationMatrix(cross,population,verbose)
  maximums <- apply(abs(cormatrix),2,max)
  means <- apply(abs(cormatrix),2,mean)
  sds <- apply(abs(cormatrix),2,sd)
  #select markers that are correlated highly with more than one of the old markers
  selected <- which(maximums > (means+corSDTreshold*sds))
  cormatrix <- cormatrix[,selected]
  bestCorMarkers <- matrix(0,length(selected),2)
  bestCorMarkers[,1] <- apply(abs(cormatrix),2,function(r){rownames(cormatrix)[which.max(r)]})
  bestCorMarkers[,2] <- apply(abs(cormatrix),2,function(r){rownames(cormatrix)[which.max(r[-which.max(r)])]})
  rownames(bestCorMarkers) <- rownames(cormatrix)
  output <- t(apply(bestCorMarkers,1,bestCorelatedSub.internal,cur_map))
  invisible(output)
}

bestCorelatedSub.internal <- function(bestCorMarkersRow,cur_map){
  chr <- cur_map[bestCorMarkersRow[1],1]
  pos <- mean(cur_map[bestCorMarkersRow[1],2],cur_map[bestCorMarkersRow[2],2])
  invisible(c(chr,pos))
}

###########################################################################################################
#                                    *** bestQTL.internal ***
#
# DESCRIPTION:
# 	subfunction of segragateChromosomes.internal, returns matrix showing for every reco map chromosome from 
#	which physicall map chromosome majority of markers comes
# OUTPUT:
#	vector with new ordering of chromosomes inside cross object
############################################################################################################
bestQTL.internal <- function(cross, population, treshold,verbose=FALSE){
  genotypes <- population$offspring$genotypes$real
  markers <- markernames(cross)
  phenotypes <- pull.geno(cross)[,markers]
  output <- NULL
  if(verbose) cat("Starting qtl analysis.\n")
  s<- proc.time()
  if(is.null(population$offspring$genotypes$qtl)) stop("No qtl data in population$offspring$genotypes$qtl, run scanQTLs function first.")
  peaksMatrix <- getpeaks.internal(abs(population$offspring$genotypes$qtl$lod),treshold)
  e<- proc.time()
  if(verbose) cat("Qtl analysis done in:",(e-s)[3],"seconds\n")
  rownames(peaksMatrix) <- markers
  for(marker in markers){
    if(sum(peaksMatrix[marker,]==2)==1){
      cur_max <- which.max(population$offspring$genotypes$qtl$lod[marker,])
      cur_row <- c(population$offspring$genotypes$qtl$chr[marker,cur_max],population$offspring$genotypes$qtl$pos[marker,cur_max])
      output <- rbind(output,cur_row)
    }else{
      output <- rbind(output,c(NA,NA))
    }
  }
  #to have same format of the output as in bestcorrelated
  rownames(output) <- markers
  if(any(is.na(output[,1]))) output <- output[-which(is.na(output[,1])),]
  if(any(is.na(output[,2]))) output <- output[-which(is.na(output[,2])),]
  invisible(output)
}


###########################################################################################################
#                                    *** QTLscan.internal ***
#
# DESCRIPTION:
# subfunction by Danny Arends to map QLTs modfied to work on a single phenotype
# OUTPUT:
#	vector with new ordering of chromosomes inside cross object
############################################################################################################
scanQTLs <- function(population,map=c("genetic","physical"),verbose=FALSE){
  if(missing(population)) stop("Please provide a population object\n")
  check.population(population)
  if(is.null(population$offspring$genotypes$real)){
    stop("No original genotypes in population$offspring$genotypes$real, load them in using intoPopulation\n")
  }
  if(is.null(population$offspring$genotypes$simulated)){
    stop("No simulated genotypes in population$offspring$genotypes$simulated, run findBiomarkers first\n")
  }
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
    returncross <- genotypesToCross.internal(population10pheno,"real","map_genetic")
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
    returncross <- genotypesToCross.internal(population10pheno,"real","map_physical")
  }
  returncross$pheno <- t(population$offspring$genotypes$simulated)
  returncross <- calc.genoprob(returncross, step=0.1)
  if(verbose) cat("Starting qtl scan, this may take a long time to finish!\n")
  s <- proc.time()
  res <- NULL
  pos <- NULL
  chr <- NULL
  names_ <- NULL
  for(i in 1:nrow(population$offspring$genotypes$simulated)){
    if(i%%50==0){
      cat("Analysing marker:",i,"\n")
      }
    curScan <- scanone(returncross,pheno.col=i,method="hk")
    pos <- rbind(pos,curScan[,2])
    res <- rbind(res,curScan[,3])
    chr <- rbind(chr,curScan[,1])
    names_ <- c(names_,colnames(returncross$pheno)[i])
  }
  #population$offspring$genotypes$qtl <- t(matrix(unlist(lapply(markers,QTLscan.internal,phenotypes,genotypes)),nrow(genotypes),length(markers)))
  e <- proc.time()
  if(verbose) cat("Qtl scan done in",(e-s)[3],"s\n")
  population$offspring$genotypes$qtl$lod <- res
  population$offspring$genotypes$qtl$pos <- pos
  population$offspring$genotypes$qtl$chr <- chr
  population$offspring$genotypes$qtl$names <- names_
  rownames(population$offspring$genotypes$qtl$lod) <- names_
  colnames(population$offspring$genotypes$qtl$lod) <- rownames(curScan)   
  rownames(population$offspring$genotypes$qtl$chr) <- names_
  colnames(population$offspring$genotypes$qtl$chr) <- rownames(curScan)  
  rownames(population$offspring$genotypes$qtl$pos) <- names_
  colnames(population$offspring$genotypes$qtl$pos) <- rownames(curScan)
  invisible(population)
}

###########################################################################################################
#                                    *** QTLscan.internal ***
#
# DESCRIPTION:
# subfunction by Danny Arends to map QLTs modfied to work on a single phenotype
# OUTPUT:
#	vector with new ordering of chromosomes inside cross object
############################################################################################################
getpeaks.internal <- function(qtlprofiles, cutoff = 4.0){
  if(!any(qtlprofiles==Inf)) qtlprofiles[which(qtlprofiles==Inf)] <- 1000
  mmatrix <- NULL
  for(x in 1:nrow(qtlprofiles)){
    peak <- FALSE
    curmax <- 0
    curmaxindex <- 1
    marker <- 1
    maximums <- NULL
    mrow <- rep(0,ncol(qtlprofiles))
    for(ab in (qtlprofiles[x,]>cutoff | qtlprofiles[x,]<(-cutoff))){
      if(ab){
        peak <- TRUE
        if(qtlprofiles[x,marker]/abs(qtlprofiles[x,marker]) > 0){
          if(qtlprofiles[x,marker] > curmax){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- marker
          }
        }else{
          if(qtlprofiles[x,marker] < (-curmax)){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- -marker
          }
        }
        if(ncol(qtlprofiles)==marker){
          if(curmax!=0) maximums <- c(maximums,curmaxindex)
        }
      }else{
        if(curmax!=0) maximums <- c(maximums,curmaxindex)
        peak <- FALSE
        curmax <- 0
      }
      marker <- marker+1
    }
    mrow[which(qtlprofiles[x,] > cutoff)] <- 1
    mrow[which(qtlprofiles[x,] < -cutoff)] <- -1
    for(a in which(maximums>0)){
      mrow[maximums[a]] <- 2
    }
    for(b in which(maximums<0)){
      mrow[(-maximums[b])] <- -2
    }
    mmatrix <- rbind(mmatrix,mrow)
  }
  mmatrix
}

###########################################################################################################
#                                    *** map2mapCorrelationMatrix.internal ***
#
# DESCRIPTION:
# 	calculating correlation matrix between genotypes inside cross object and ones from population
# OUTPUT:
#   matrix of correlations
############################################################################################################
map2mapCorrelationMatrix<- function(cross,population,verbose=FALSE){
  if(missing(cross)) stop("Please provide a cross object\n")
  if(missing(population)) stop("Please provide original genotypes\n")
  #is.cross(cross)
  genotypes <- pull.geno(cross)
  if(verbose) cat("Calculating correlation matrix\n")
  if(!is.null(population$offspring$genotypes$real)){
    genotypesCorelationMatrix <- apply(genotypes,2,function(cgc){cor(cgc,t(population$offspring$genotypes$real),use="pair")})
    colnames(genotypesCorelationMatrix) <- colnames(genotypes)
    rownames(genotypesCorelationMatrix) <- rownames(population$offspring$genotypes$real)    
    invisible(genotypesCorelationMatrix)
  }else{
    stop("Load known genotypes into the population using intoPopulation(p,genotypes,\"offspring$genotypes\")")
  }
}

###########################################################################################################
#                                           *** map2mapImage ***
#
# DESCRIPTION:
# 	subfunction of segragateChromosomes.internal, returns matrix showing for every reco map chromosome from 
#	which physicall map chromosome majority of markers comes# 
# OUTPUT:
#	vector with new ordering of chromosomes inside cross object
############################################################################################################
map2mapImage <- function(genotypesCorelationMatrix,population,cross,corThreshold=0.5,verbose=FALSE){
  if(missing(genotypesCorelationMatrix)){
    cat("Correlation matrix not provided, calulating one")
    genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross,population,verbose)
    if(missing(cross)) stop("No object of class cross, please run either cross.denovo or enrichExistingMap\n")
    if(missing(population)) stop("Please provide a population object\n")
    check.population(population)
  }
  heatmap(genotypesCorelationMatrix,breaks = c(-1,-corThreshold,corThreshold,1),col=c("blue","white","red"),Rowv=NA)
}
