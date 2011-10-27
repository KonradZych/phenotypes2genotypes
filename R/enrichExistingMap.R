############################################################################################################
#
# rearrangeMarkers.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified October 2011
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
# Contains: enrichExistingMap, rearrangeMarkers,
#             bestCorelated.internal, map2mapCorrelationMatrix.internal, map2mapImage
#
############################################################################################################

###########################################################################################################
#                                    *** enrichExistingMap ***
#
# DESCRIPTION:
# 	enriching existing genetic map adding markers derived from gene expression
# 
# PARAMETERS:
# 	cross - object of class cross, containing physical or genetic map
# 	map - which map should be used for comparison:
#			- genetic - genetic map from cross$maps$genetic
#			- physical - physical map from cross$maps$physical
#	corTreshold - markers not having corelation above this number with any of chromosomes are removed
#	verbose - be verbose
#
#
# OUTPUT:
#	object of class cross
#
############################################################################################################
enrichExistingMap <- function(population,cross,map=c("genetic","physical"),corTreshold=0.6,verbose=FALSE,debugMode=0){
  if(missing(population)) stop("Please provide a population object\n")
  if(is.null(population$offspring$genotypes$real)){
    stop("No original genotypes in population$offspring$genotypes$real, load them in using intoPopulation\n")
  }
  check.population(population)
  
  if(missing(cross)){
    if(is.null(population$offspring$genotypes$simulated)){
      stop("No genotype data in population$offspring$genotypes$simulated, run toGenotypes first\n")
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
	cross <- rearrangeMarkers(cross,population,map,corTreshold,TRUE,verbose=verbose)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Enrichment of original map done in:",(e1-s1)[3],"seconds.\n")
  
  #*******ORDERING NEW MAP*******
  if(verbose)cat("Ordering markers inside the cross object\n")
  s1 <- proc.time()
  aa <- tempfile()
  sink(aa)
  cross <- orderMarkers(cross,use.ripple=F,verb=T)
  sink()
  file.remove(aa)
  e1 <- proc.time()
  if(verbose && debugMode==2)cat("Saving data into cross object done in:",(e1-s1)[3],"seconds.\n")
	invisible(cross)
}


###########################################################################################################
#                                    *** rearrangeMarkers ***
#
# DESCRIPTION:
# 	ordering chromosomes using genetic/physical map and corelation rule
# 
# PARAMETERS:
# 	cross - object of class cross, containing physical or genetic map
# 	map - which map should be used for comparison:
#			- genetic - genetic map from cross$maps$genetic
#			- physical - physical map from cross$maps$physical
#	corTreshold - markers not having corelation above this number with any of chromosomes are removed
#	addMarkers - should markers used for comparison be added to output cross object
#	verbose - be verbose
#
#
# OUTPUT:
#	object of class cross
#
############################################################################################################
rearrangeMarkers <- function(cross,population,map=c("genetic","physical"),corTreshold=0.6,addMarkers=FALSE,verbose=FALSE){
  if(missing(cross)) stop("Please provide a cross object\n")
  if(missing(population)) stop("Please provide a population object\n")
  check.population(population)
  output <- bestCorelated.internal(cross,population,corTreshold,verbose)
  if(verbose) cat("selected",nrow(output),"markers for further analysis\n")
  map <- defaultCheck.internal(map,"map",2,"genetic")
	if(map=="genetic"){
    cur_map <- population$maps$genetic
  }else{
    cur_map <- population$maps$physical
  }
  if(verbose) cat("old map contains",max(cur_map[,1]),"chromosomes\n")
	cross_ <- cross
	cross_$geno <- vector(max(cur_map[,1]), mode="list")
	cross_$pheno <- pull.pheno(cross)
	if(verbose) cat("Reordering markers \n")  
	for(x in 1:max(cur_map[,1])){
		if(verbose) cat("- chr ",x," -\n")    
		oldnames <- rownames(cur_map)[which(cur_map[,1]==x)]
    newmarkers <- which(output[,2]%in%oldnames)
    if(verbose) cat("Selected:",length(newmarkers),"new and",length(oldnames),"original markers \n") 
		if(addMarkers){
			cross_$geno[[x]]$data <- cbind(pull.geno(cross)[,output[newmarkers,1]],t(population$offspring$genotypes$real[oldnames,]))
			newmap <- 1:(length(newmarkers)+length(oldnames))
			names(newmap) <- c(output[newmarkers,1],oldnames)
		}else{
			cross_$geno[[x]]$data <- pull.geno(cross)[,output[newmarkers,1]]
			newmap <- 1:length(newmarkers)
			names(newmap) <- output[newmarkers,1]
		}
		cross_$geno[[x]]$map <- c(newmap)
	}
	names(cross_$geno) <- 1:length(cross_$geno)
	for(i in 1:length(cross_$geno)){
		class(cross_$geno[[i]]) <- "A"
	}
	invisible(cross_)
}


###########################################################################################################
#                                    *** bestCorelated.internal ***
#
# DESCRIPTION:
# 	subfunction of segragateChromosomes.internal, returns matrix showing for every reco map chromosome from 
#	which physicall map chromosome majority of markers comes
# 
# PARAMETERS:
# 	cross - object of class cross, containing physical or genetic map
#	cur_map - object containing map to be used
#	verbose - be verbose
# 
# OUTPUT:
#	vector with new ordering of chromosomes inside cross object
#
############################################################################################################
bestCorelated.internal <- function(cross,population,corTreshold,verbose=FALSE){
  genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross,population,verbose)
  #select markers that are correlated highly with more than one of the old markers
  selected <- which(apply(abs(genotypesCorelationMatrix),2,function(r){length(which(r > corTreshold))})!=0)
  genotypesCorelationMatrix <- genotypesCorelationMatrix[,selected]
  output <- matrix(0,length(selected),2)
  output[,1] <- colnames(genotypesCorelationMatrix)
  output[,2] <- apply(abs(genotypesCorelationMatrix),2,function(r){rownames(genotypesCorelationMatrix)[which.max(r)]})
  rownames(output) <- colnames(genotypesCorelationMatrix)
  invisible(output)
}

###########################################################################################################
#                                    *** map2mapCorrelationMatrix.internal ***
#
# DESCRIPTION:
# 	calculating correlation matrix between genotypes inside cross object and ones from population
# 
# PARAMETERS:
#   cross - an object of class cross
#   population - an object of class population
#   verbose - be verbose
# 
# OUTPUT:
#   matrix of correlations
#
############################################################################################################
map2mapCorrelationMatrix<- function(cross,population,verbose=FALSE){
  if(missing(cross)) stop("Please provide a cross object\n")
  if(missing(population)) stop("Please provide a population object\n")
  check.population(population)
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
#	which physicall map chromosome majority of markers comes
# 
# PARAMETERS:
# 	cross - object of class cross, containing physical or genetic map
#	cur_map - object containing map to be used
#	verbose - be verbose
# 
# OUTPUT:
#	vector with new ordering of chromosomes inside cross object
#
############################################################################################################
map2mapImage <- function(genotypesCorelationMatrix,population,cross,corThreshold=0.5,verbose=FALSE){
  if(missing(genotypesCorelationMatrix)){
    cat("Correlation matrix not provided, calulating one")
    genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross,population,verbose)
    if(missing(cross)) stop("No object of class cross, please run either createNewMap or enrichExistingMap\n")
    if(missing(population)) stop("Please provide a population object\n")
    check.population(population)
  }
  heatmap(genotypesCorelationMatrix,breaks = c(-1,-corThreshold,corThreshold,1),col=c("blue","white","red"),Rowv=NA)
}
