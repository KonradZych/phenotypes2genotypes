############################################################################################################
#
# toGenotypes.R
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
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
#
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: toGenotypes
#           convertToGenotypes.internal, splitPheno.internal, selectMarkersUsingMap.internal,
#			filterGenotypes.internal, filterRow.internal, splitRowSubEM.internal
#
############################################################################################################

############################################################################################################
#									*** toGenotypes ***
#
# DESCRIPTION:
#	function that chooses from the matrix only appropriate markers with specified rules
# 
# PARAMETERS:
# 	population - Ril type object, must contain founders phenotypic data.
# 	orderUsing- which map should be used to order markers (default - none)
# 		- map_genetic - genetic map
#		- map_physical - physical map
# 	treshold - If Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.
# 	overlapInd - Number of individuals that are allowed in the overlap
# 	proportion - Proportion of individuals expected to carrying a certain genotype 
# 	margin - Proportion is allowed to varry between this margin (2 sided)
# 	minChrLength -if maximal distance between the markers in the chromosome is lower than this value,
#		whole chromosome will be dropped
# 	verbose - Be verbose
# 	debugMode - 1: Print our checks, 2: print additional time information
# 
# OUTPUT:
#	an object of class cross
#
############################################################################################################
toGenotypes <- function(population, treshold=0.05, overlapInd = 0, proportion = c(50,50), margin = 15, verbose=FALSE, debugMode=0){
	#*******CHECKS*******
	check.population(population)
	s<-proc.time()
	if(any(proportion < 1) || sum(proportion) != 100) stop("Wrong proportion paramete\n")
	if(any(!(is.numeric(population$founders$phenotypes)))){
		population <- intoPopulation(population, population$founders$phenotypes, "founders")
	}
	if(any(!(is.numeric(population$offspring$phenotypes)))){
		population <- intoPopulation(population, population$offspring$phenotypes, "offspring$phenotypes")
	}
	if(overlapInd < 0 || overlapInd > ncol(population$offspring$phenotypes)) stop("overlapInd is a number (0,lenght of the row).")
	if(verbose && debugMode==1) cat("toGenotypes starting withour errors in checkpoint.\n")
	
	#*******CONVERTING CHILDREN PHENOTYPIC DATA TO GENOTYPES*******
	s1 <- proc.time()
	population <- convertToGenotypes.internal(population, treshold, overlapInd, proportion, margin, verbose, debugMode)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Converting phenotypes to genotypes done in:",(e1-s1)[3],"seconds.\n")
	
	#*******RETURNING CROSS OBJECT*******
	e<-proc.time()
	if(verbose) cat("toGenotypes done in",(e-s)[3],"seconds\n")
	invisible(population)
}

############################################################################################################
#									*** convertToGenotypes.internal ***
#
# DESCRIPTION:
#	function splitting differentially expressed markers into two genotypes
# 
# PARAMETERS:
# 	population - object of class population, must contain founders phenotypic data.
# 	orderUsing- which map should be used to order markers (default - none)
# 		- map_genetic - genetic map
#		- map_physical - physical map
# 	treshold - if Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.
# 	overlapInd - number of individuals that are allowed in the overlap
# 	proportion - proportion of individuals expected to carrying a certain genotype 
# 	margin - proportion is allowed to varry between this margin (2 sided)
# 	verbose - be verbose
# 	debugMode - 1: Print our checks, 2: print additional time information 
# 
# OUTPUT:
#	object of class population
#
############################################################################################################
convertToGenotypes.internal <- function(population, treshold, overlapInd, proportion, margin, verbose=FALSE, debugMode=0){
	### initialization
	if(verbose && debugMode==1) cat("convertToGenotypes starting.\n")
	output <- NULL
	markerNames <- NULL 
	
	### selection step
	### up-regulated
	upNotNull <- which(population$founders$RP$pval[,1] > 0)
	upBelowTreshold <- which(population$founders$RP$pval[,1] < treshold)
	upSelected <- upBelowTreshold[which(upBelowTreshold%in%upNotNull)]
	upParental <- population$founders$phenotypes[upSelected,]
	upRils <- population$offspring$phenotypes[rownames(upParental),]
	### down-regulated
	downNotNull <- which(population$founders$RP$pval[,2] > 0)
	downBelowTreshold <- which(population$founders$RP$pval[,2] < treshold)
	downSelected <- downBelowTreshold[which(downBelowTreshold%in%downNotNull)]
	downParental <- population$founders$phenotypes[downSelected,]
	downRils <- population$offspring$phenotypes[rownames(downParental),]
	
	### checking if anything is selected and if yes - processing
	if(!(is.null(dim(upRils)))&&(nrow(upRils)!=0)){
		if(!(is.null(dim(downRils)))&&(nrow(downRils)!=0)){
			# best situation
			if(verbose) cat("Selected",nrow(downRils),"markers (DOWN), ",nrow(upRils),"markers (UP).\n")
			inupndown <- which(rownames(upRils) %in% rownames(downRils))
			if(verbose&&length(inupndown)>0){
        cat("WARNING: Overlap between UP n DOWN:",length(inupndown),", removing from UP.\n")
        upRils <- upRils[-inupndown,]
      }
			cur <- splitPheno.internal(downRils, downParental, overlapInd, proportion, margin, population$founders$groups, 0, nrow(upRils),verbose)
			output <- rbind(output,cur[[1]])
			markerNames <- c(markerNames,cur[[2]])
		}else{
			if(verbose) cat("Selected ",nrow(upRils),"upregulated markers.\n")
		}
		cur <- splitPheno.internal(upRils, upParental, overlapInd, proportion, margin, population$founders$groups, 1, 0,verbose)
		output <- rbind(output,cur[[1]])
		markerNames <- c(markerNames,cur[[2]])
	}else{
		if(!(is.null(dim(downRils)))&&(nrow(downRils)!=0)){
			if(verbose) cat("Selected ",nrow(downRils),"downregulated markers.\n")
			cur <- splitPheno.internal(downRils, downParental, overlapInd, proportion, margin, population$founders$groups, 0,0,verbose)
			output <- rbind(output,cur[[1]])
			markerNames <- c(markerNames,cur[[2]])
		}else{
			stop("None of the markers was selected using specified treshold: ",treshold,"\n")
		}
	}
	
	### putting results inside population object
	if(is.null(dim(output))) stop("No markers selected.")
	population$offspring$genotypes$simulated <- output
	colnames(population$offspring$genotypes$simulated) <- colnames(upRils)
	rownames(population$offspring$genotypes$simulated) <- markerNames
	invisible(population)
}

############################################################################################################
#									*** splitPheno.internal ***
#
# DESCRIPTION:
#	subfunction of convertToGenotypes.internal, splitting children markers using founders mean values
# 
# PARAMETERS:
# 	offspring - matrix of up/down regulated genes in offspring
# 	founders - matrix of up/down regulated genes in parents
# 	overlapInd - Number of individuals that are allowed in the overlap
# 	proportion - Proportion of individuals expected to carrying a certain genotype 
# 	margin - Proportion is allowed to varry between this margin (2 sided)
# 	groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
# 	up - 1 - genes up 0 - down regulated
# 
# OUTPUT:
#	list containg genotype matrix and names of selected markers
#
############################################################################################################
#DANNY: TODO MERGE splitPhenoRowEM.internal into this function 
##K: left, I\'ll try to use apply here instead of for
splitPheno.internal <- function(offspring, founders, overlapInd, proportion, margin, groupLabels, up, left=0, verbose=FALSE){
	output <- NULL
	markerNames <- NULL
  s <-proc.time()
	for(x in 1:nrow(offspring)){
		cur <- splitPhenoRowEM.internal(x, offspring, founders, overlapInd, proportion, margin, groupLabels, up, verbose)
		if(!(is.null(cur))){
			output <- rbind(output,cur)
			markerNames <- c(markerNames,rownames(offspring)[x])
		}
    if(x%%100==0){
      e <- proc.time()
      te <- ((e-s)[3]/x)*(nrow(offspring)-x+left)
      cat("Done with marker",x,"/",nrow(offspring)+left,". Time remaining:",te,"s\n")
    }
	}
	invisible(list(output,markerNames))
}

############################################################################################################
#									*** splitPhenoRowEM.internal ***
#
# DESCRIPTION:
#	subfunction of splitRow.internal, splitting one row using EM algorithm
# 
# PARAMETERS:
# 	x - name of currently processed row
# 	offspring - matrix of up/down regulated genes in offspring
# 	founders - matrix of up/down regulated genes in parents
# 	overlapInd - Number of individuals that are allowed in the overlap
# 	proportion - Proportion of individuals expected to carrying a certain genotype 
# 	margin - Proportion is allowed to varry between this margin (2 sided)
# 	groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
# 	up - 1 - genes up 0 - down regulated
# 
# OUTPUT:
#	genotype row
#
############################################################################################################
splitPhenoRowEM.internal <- function(x, offspring, founders, overlapInd, proportion, margin, groupLabels, up=1,verbose=FALSE){
	y<-x
	x<-rownames(offspring)[x]
	aa <- tempfile()
	sink(aa)
	nrDistributions <- length(proportion)
	result <- rep(0,length(offspring[x,]))
	
	EM <- NULL
	s1<-proc.time()
	tryCatch(EM <- normalmixEM(sort(offspring[x,]), k=nrDistributions, maxrestarts=0, maxit = 100,fast=FALSE),error = function(x){cat(x[[1]],"\n")})
	e1<-proc.time()

	if(is.null(EM)){
	 result <- NULL
	}else{
		if(up==1){
			genotypes <- c(0:(nrDistributions-1))
		}else if(up==0){
			genotypes <- c((nrDistributions-1):0)
		}
		
		len <- vector(mode="numeric",length=nrDistributions)
		 for(i in 1:nrDistributions){
			len[i]<-length(offspring[x,])*EM$lambda[which(EM$mu==sort(EM$mu)[i])]
			startVal <- sum(len[1:i-1])
			result[which(offspring[x,] %in% sort(offspring[x,])[startVal:(startVal+len[i])])] <- genotypes[i]
		 }
     #Danny: WHAT IS THIS ????, I removed the 2 functions
		 #if(!checkMu.internal(offspring,EM,overlapInd)){
     #  result <- middleDistribution.internal(offspring,result,EM)
		 #}
      result <- filterRow.internal(result, overlapInd, proportion, margin, genotypes)
	}
	sink()
	file.remove(aa)
	invisible(result)
}

############################################################################################################
#									*** filterRow.internal ***
#
# DESCRIPTION:
#	checking whether given genotype row meets given requirments
# 
# PARAMETERS:
# 	population - Ril type object, must contain founders phenotypic data.
# 	overlapInd - Number of individuals that are allowed in the overlap
# 	proportion - Proportion of individuals expected to carrying a certain genotype 
# 	margin - Proportion is allowed to varry between this margin (2 sided)
# 	verbose - Be verbose
# 	debugMode - 1: Print our checks, 2: print additional time information
# 
# OUTPUT:
#	genotype row if it meets requirments or NULL
#
############################################################################################################
filterRow.internal <- function(result, overlapInd, proportion, margin, genotypes){
	### creating inverted genotypes matrix, to be sure, that we won't filter out anythng in correct proportion
	### this function returns either unchanged result vector, which is then rbinded to other results, or
	### NULL, which is ignored by rbind, and we drop current result
	genotypes2 <- genotypes[length(genotypes):1]
	if(filterRowSub.internal(result, overlapInd, proportion, margin, genotypes)){
		invisible(result)
	}else if(filterRowSub.internal(result, overlapInd, proportion, margin, genotypes2)){
		invisible(result)
	}else{
		return(NULL)
	}
}

############################################################################################################
#									*** filterRowSub.internal ***
#
# DESCRIPTION:
# 	subfunction of filterGenotypes.internal, filtering one row
# 
# PARAMETERS:
# 	genotypeRow - currently processed row
# 	overlapInd - Number of individuals that are allowed in the overlap
# 	proportion - Proportion of individuals expected to carrying a certain genotype 
# 	margin - Proportion is allowed to varry between this margin (2 sided)
# 
# OUTPUT:
#	boolean
#
############################################################################################################
filterRowSub.internal <- function(genotypeRow, overlapInd, proportion, margin, genotypes){
	if(sum(is.na(genotypeRow)) > overlapInd) return(FALSE)
	for(i in 1:length(proportion)){
		cur <- sum(genotypeRow==genotypes[i])/length(genotypeRow) * 100
		if(is.na(cur)) return(FALSE)
		upLimit <- proportion[i] + margin/2
		downLimit <- proportion[i] - margin/2
		if(!((cur < upLimit)&&(cur > downLimit))){
			return(FALSE)
		}
	}
	return(TRUE)
}
