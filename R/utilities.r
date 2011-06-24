############################################################################################################
#
# utilities.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified June 2011
# last modified in version: 0.7.2
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
# Contains: toGenotypes
#           convertToGenotypes.internal, splitRow.internal, filterGenotypes.internal, filterRow.internal
#			sortMap.internal, intoRil, removeChromosomes.internal
#
############################################################################################################

############################################################################################################
#cleanMap - removing markers that cause reco map to be too long
# 
# cross - R/qtl cross type object
# difPercentage - If removing marker will make map shorter by this percentage of its length, then it will be dropped.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information 
#
############################################################################################################
cleanMap <- function(cross, difPercentage, verbose=FALSE, debugMode=0){
	if(verbose && debugMode==1) cat("cleanMap starting withour errors in checkpoint.\n")
	s <- proc.time()
	for(i in 1:length(cross$geno)){
		begMarkers <- length(cross$geno[[i]]$map)
		begLength <- max(cross$geno[[i]]$map)
		for(j in names(cross$geno[[i]]$map)){
			if(max(cross$geno[[i]]$map)>150){
				cur_max <- max(cross$geno[[i]]$map)
				cross2 <- drop.markers(cross,j)
				newmap <- est.map(cross2,offset=0)
				cross2 <- replace.map(cross2, newmap)
				new_max <- max(cross2$geno[[i]]$map)
				dif <- cur_max-new_max
				if(dif > (difPercentage/100 * cur_max)){
					if(verbose) cat("------Removed marker:",j,"to make chromosome",i,"map smaller from",cur_max,"to",new_max,"\n")
					cross <- cross2
				}
			}
		}
		removed <- begMarkers-length(cross$geno[[i]]$map)
		if(removed>0)cat("Removed",removed,"out of",begMarkers,"markers on chromosome",i," which led to shortening map from ",begLength,"to",max(cross$geno[[i]]$map),"(",100*(begLength-max(cross$geno[[i]]$map))/begLength,"%)\n")
	}
	e <- proc.time()
	if(verbose && debugMode==2)cat("Map cleaning done in:",(e-s)[3],"seconds.\n")
	invisible(cross)
}

############################################################################################################
#print.population - overwritting print function for objects of class "population"
# 
# x - object of class population
# ... - passed to cats
#
############################################################################################################
print.population <- function(x,...){
	cat("*************************************************************************************\n")
	cat("This is object of class population, too complex to print it, so we provide you with summary.\n")
	if(!(is.null(x$offspring))){
		if(!(is.null(x$offspring$phenotypes))){
			cat("Object contains phenotypic data for",ncol(x$offspring$phenotypes),"offspring individuals covering",nrow(x$offspring$phenotypes),"probes.\n",...)
		}else{
			stop("There is no phenotypic data for offspring in this object. This is not acceptable in real ril object.\n",...)
		}
		if(!(is.null(x$offspring$genotypes$read))){
			cat("Object contains genotypic data for",ncol(x$offspring$phenotypes),"offspring individuals covering",nrow(x$offspring$phenotypes),"probes.\n",...)
		}else{
			cat("There is no genotypic data for offspring in this object.\n",...)
		}
		if(!(is.null(x$offspring$map))){
			cat("Object contains physical map covering",nrow(x$offspring$map),"markers from",length(table(x$offspring$map[,1])),"chromosomes.\n",...)
		}else{
			cat("There is no physical genetic map in this object.\n")
		}
	}else{
		cat("WARNING: There is no phenotypic data for offspring. This is not acceptable in real ril object.\n",...)
	}
	
	if(!(is.null(x$founders))){
		if(!(is.null(x$founders$phenotypes))){
			cat("Object contains phenotypic data for",ncol(x$founders$phenotypes),"founders individuals covering",nrow(x$founders$phenotypes),"probes.\n",...)
		}else{
			stop("There is no phenotypic data for parents in this object. This is not acceptable in real ril object.\n",...)
		}
		if(!(is.null(x$founders$RP))){
			cat("Object contains RP analysis results.\n",...)
		}else{
			cat("There is no RP analysis result in this object.\n",...)
		}
		if(!(is.null(x$founders$groups))){
			cat("Parental groups are as following:",x$founders$groups,"\n",...)
		}else{
			cat("There is no information about founders groups in this object.\n",...)
		}
	}else{
		cat("WARNING: There is no phenotypic data for parents. This is not acceptable in real ril object.\n",...)
	}
	cat("*************************************************************************************\n")
}

############################################################################################################
#removeChromosomes.internal: removing from cross chromosomes that have too little markers
# 
# cross - object of R/qtl cross type
# minChrLength - minimal number of markers chromosome have to contaion not to be removed)
#
############################################################################################################
removeChromosomes.internal <- function(cross, numberOfChromosomes){
	for(i in length(cross$geno):(numberOfChromosomes+1)){
		cat("removing chromosome:",i," markers:",names(cross$geno[[i]]$map),"\n")
		cross$rmv <- cbind(cross$rmv,cross$geno[[i]]$data)
		cross <- drop.markers(cross, names(cross$geno[[i]]$map))
		names(cross$geno) <- 1:length(cross$geno)
	}
	invisible(cross)
}
