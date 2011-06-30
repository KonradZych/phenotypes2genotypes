############################################################################################################
#
# toGenotypes.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified June 2011
# last modified in version: 0.8.1 
# in current version: active, in main workflow
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
#           sortMap.internal, majorityRule.internal, mergeChromosomes.internal, 
#
############################################################################################################


############################################################################################################
#									*** segragateChromosomes.internal ***
#
# DESCRIPTION:
# 	ordering chromosomes using physical map and majority rule
# 
# PARAMETERS:
# 	cross - object of class cross, containing physical or genetic map
# 
# OUTPUT:
#	object of class cross
#
############################################################################################################
segregateChromosomes.internal <- function(cross){
	if(is.null(cross$maps$physical)){
		cat("WARNING: no physical map, function will return unchanged cross object\n")
	}else{
		output <- majorityRule.internal(cross)
		print(output)
		### until every chr on phys map is match exactly once
		while(max(apply(output,2,sum))>1){
			toMerge <- which(apply(output,2,sum)>1)
			for(curToMerge in toMerge){
				curMerge <- which(output[,curToMerge]==max(output[,curToMerge]))
				map <- cross$maps$physical
				cross <- mergeChromosomes.internal(cross,curMerge,curMerge[1])
				cross$maps$physical <- map
				output <- majorityRule.internal(cross)
			}
		}
		
		order1 <- matrix(0,ncol(output),nrow(output))
		order2 <- matrix(1,ncol(output),nrow(output))
		### until next iteration doesn't change the result
		while(any(order1!=order2)){
			order1 <- output
			for(l in 1:ncol(output)){
				cur <- which(output[,l]==max(output[,l]))
				if(cur!=l)cross <- switchChromosomes.internal(cross,cur,l)
				output <- majorityRule.internal(cross)
			}
			order2 <- output
		}
		names(cross$geno) <- 1:length(cross$geno)
	}
	invisible(cross)
}

############################################################################################################
#									*** majorityRule.internal ***
#
# DESCRIPTION:
# 	subfunction of segragateChromosomes.internal, returns matrix showing for every reco map chromosome from 
#	which physicall map chromosome majority of markers comes
# 
# PARAMETERS:
# 	cross - object of class cross, containing physical or genetic map
# 
# OUTPUT:
#	vector with new ordering of chromosomes inside cross object
#
############################################################################################################
majorityRule.internal <- function(cross){
	knchrom <- length(table(cross$maps$physical[[1]][,1]))
	result <- matrix(0, length(cross$geno), knchrom)
	output <- matrix(0, length(cross$geno), knchrom)
	for(i in 1:length(cross$geno)){
		cur_ys <- colnames(cross$geno[[i]]$data)
		cur_xs <- cross$maps$physical[[1]][cur_ys,]
		for(j in 1:knchrom){
			result[i,j] <- sum(cur_xs[,1]==j)/nrow(cur_xs)
		}
		output[i,which(result[i,]==max(result[i,]))] <- 1
	}
	rownames(result) <- 1:nrow(result)
	colnames(result) <- 1:ncol(result)
	rownames(output) <- 1:nrow(output)
	colnames(output) <- 1:ncol(output)
	
	if(min(apply(output,2,max))!=1){
		toCheck <- which(apply(output,2,sum)!=1)
		for(x in toCheck){
			output[,x] <- 0
			output[which(result[,x]==max(result[,x])),x] <- 1
		}
	}	
	invisible(output)
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
# 	name - chromosome
# 
# OUTPUT:
#	object of class cross
#
############################################################################################################
mergeChromosomes.internal <- function(cross, chromosomes, name){
	cat("Merging chromosomes",chromosomes,"to form chromosome",name,"names:",names(cross$geno),"\n")
	geno <- cross$geno
	markerNames <- NULL
	for(j in chromosomes){
		if(j!=name) markerNames <- c(markerNames, colnames(geno[[j]]$data))
	}
	for(k in markerNames) cross <- movemarker(cross, k, name)
	cat("Ordering markers on newly merged chromosome\n")
	#cross <- orderMarkers(cross, chr=name)
	invisible(cross)
}