############################################################################################################
#
# toGenotypes.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified July 2011
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
#           convertToGenotypes.internal, splitPheno.internal, selectMarkersUsingMap.internal,
#			filterGenotypes.internal, filterRow.internal, splitRowSubEM.internal, checkMu.internal
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
# 	genotype - Which genotypic matrix should be saved to file, 
#		- simulated - made by toGenotypes
# 		- real - supported by user and read from file, 
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
#	object of class cross
#
############################################################################################################
toGenotypes <- function(population, genotype=c("simulated","real"), orderUsing=c("none","map_genetic","map_physical"),treshold=0.01, overlapInd = 0, proportion = c(50,50), margin = 15, verbose=FALSE, debugMode=0){
	#*******CHECKS*******
	is.population(population)
	s<-proc.time()
	if(proportion < 1 || proportion > 99) stop("Proportion is a percentage (1,99)")
	if(any(!(is.numeric(population$founders$phenotypes)))){
		population <- intoPopulation(population, population$founders$phenotypes, "founders")
	}
	if(any(!(is.numeric(population$offspring$phenotypes)))){
		population <- intoPopulation(population, population$offspring$phenotypes, "offspring$phenotypes")
	}
	if(overlapInd < 0 || overlapInd > ncol(population$offspring$phenotypes)) stop("overlapInd is a number (0,lenght of the row).")
	if(margin < 0 || margin > proportion) stop("Margin is a percentage (0,proportion)")
	if(verbose && debugMode==1) cat("toGenotypes starting withour errors in checkpoint.\n")
	inListCheck.internal(genotype,"genotype",c("simulated","real"))
	inListCheck.internal(orderUsing,"orderUsing",c("none","map_genetic","map_physical"))
	
	
	#*******CONVERTING CHILDREN PHENOTYPIC DATA TO GENOTYPES*******
	if(genotype=="simulated"){
		s1 <- proc.time()
		population <- convertToGenotypes.internal(population, orderUsing, treshold, overlapInd, proportion, margin, verbose, debugMode)
		e1 <- proc.time()
		if(verbose && debugMode==2)cat("Converting phenotypes to genotypes done in:",(e1-s1)[3],"seconds.\n")
	}
	
	
	#*******SAVING CROSS OBJECT*******
	s1 <- proc.time()
	cross <- genotypesToCross.internal(population,genotype=genotype,orderUsing=orderUsing,verbose=verbose,debugMode=debugMode)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Creating cross object done in:",(e1-s1)[3],"seconds.\n")
	
	#*******ADDING MAPS TO THE CROSS*******
	if(!(is.null(population$maps$physical))) cross$maps$physical <- population$maps$physical
	if(!(is.null(population$maps$genetic))) cross$maps$genetic <- population$maps$genetic
	
	#*******ADDING REAL GENOTYPE TO THE CROSS*******
	if(!(is.null(population$offspring$genotypes$real))&&genotype=="simulated") cross$genotypes$real <- population$offspring$genotypes$real
	
	#*******RETURNING CROSS OBJECT*******
	e<-proc.time()
	if(verbose) cat("toGenotypes done in",(e-s)[3],"seconds\n")
	invisible(cross)
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
convertToGenotypes.internal <- function(population, orderUsing, treshold, overlapInd, proportion, margin, verbose=FALSE, debugMode=0){
	### initialization
	if(verbose && debugMode==1) cat("convertToGenotypes starting.\n")
	output <- NULL
	markerNames <- NULL 
	
	### selection step
	### up-regulated
	upNotNull <- which(population$founders$RP$pval[1] > 0)
	upBelowTreshold <- which(population$founders$RP$pval[1] < treshold)
	upSelected <- upBelowTreshold[which(upBelowTreshold%in%upNotNull)]
	upParental <- population$founders$phenotypes[upSelected,]
	upParental <- selectMarkersUsingMap.internal(upParental,population,orderUsing,verbose,debugMode)
	upRils <- population$offspring$phenotypes[rownames(upParental),]
	### down-regulated
	downNotNull <- which(population$founders$RP$pval[2] > 0)
	downBelowTreshold <- which(population$founders$RP$pval[2] < treshold)
	downSelected <- downBelowTreshold[which(downBelowTreshold%in%downNotNull)]
	downParental <- population$founders$phenotypes[downSelected,]
	downParental <- selectMarkersUsingMap.internal(downParental,population,orderUsing,verbose,debugMode)
	downRils <- population$offspring$phenotypes[rownames(downParental),]
	
	### checking if anything is selected and if yes - processing
	if(!(is.null(dim(upRils)))){
		if(!(is.null(dim(downRils)))){
			# best situation
			if(verbose) cat("Selected ",nrow(upRils),"upregulated markers and ",nrow(downRils),"downregulated markers.\n")
			cur <- splitPheno.internal(downRils, downParental, overlapInd, proportion, margin, population$founders$groups, 0)
			output <- rbind(output,cur[[1]])
			markerNames <- c(markerNames,cur[[2]])
		}else{
			if(verbose) cat("Selected ",nrow(upRils),"upregulated markers.\n")
		}
		cur <- splitPheno.internal(upRils, upParental, overlapInd, proportion, margin, population$founders$groups, 1)
		output <- rbind(output,cur[[1]])
		markerNames <- c(markerNames,cur[[2]])
	}else{
		if(!(is.null(dim(downRils)))){
			if(verbose) cat("Selected ",nrow(downRils),"downregulated markers.\n")
			cur <- splitPheno.internal(downRils, downParental, overlapInd, proportion, margin, population$founders$groups, 0)
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
#									*** selectMarkersUsingMap.internal ***
#
# DESCRIPTION:
#	selecting from phenpotypeMatrix only markers present on map selected for ordering
# 
# PARAMETERS:
# 	phenotypeMatrix - matrix - rows - markers, cols - individuals
# 	population - object of class population, must contain founders phenotypic data.
# 	orderUsing- which map should be used to order markers (default - none)
# 		- map_genetic - genetic map
#		- map_physical - physical map
# 	verbose - be verbose
# 	debugMode - 1: Print our checks, 2: print additional time information 
# 
# OUTPUT:
#	object of class population
#
############################################################################################################
selectMarkersUsingMap.internal <- function(phenotypeMatrix,population,orderUsing,verbose,debugMode){
	if(orderUsing!="none"){
		if(orderUsing=="map_genetic"){
			if(any(!(rownames(phenotypeMatrix)%in%rownames(population$maps$genetic)))){
				cat("Not all markers selected for genotype reconstruction are present on genetic map.\n")
				phenotypeMatrix <- mapMarkers.internal(phenotypeMatrix,population$maps$genetic,1,verbose,debugMode)
			}
		}else if(orderUsing=="map_physical"){
			if(any(!(rownames(phenotypeMatrix)%in%rownames(population$maps$physical)))){
				cat("Not all markers selected for genotype reconstruction are present on physical map.\n")
				phenotypeMatrix <- mapMarkers.internal(phenotypeMatrix,population$maps$physical,1,verbose,debugMode)
			}
		}
	}
	invisible(phenotypeMatrix)
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
splitPheno.internal <- function(offspring, founders, overlapInd, proportion, margin, groupLabels, up){
	output <- NULL
	markerNames <- NULL
	for(x in rownames(offspring)){
		cur <- splitPhenoRowEM.internal(x, offspring, founders, overlapInd, proportion, margin, groupLabels, up)
		if(!(is.null(cur))){
			output <- rbind(output,cur)
			markerNames <- c(markerNames,x)
		}
	}
	invisible(list(output,markerNames))
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
	if(length(proportion)==2){
		genotypes2 <- genotypes[c(2,1)]
	}else if(length(proportion)==3){
		genotypes2 <- genotypes[c(3,2,1)]
	}
	if(filterRowSub.internal(result, overlapInd, proportion, margin, genotypes)){
		invisible(result)
	}else if(filterRowSub.internal(result, overlapInd, proportion, margin, genotypes)){
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
splitPhenoRowEM.internal <- function(x, offspring, founders, overlapInd, proportion, margin, groupLabels, up=1){
	### initialization
	print(x)
	downLimit <- mean(offspring[x,]) - 2*sd(offspring[x,])
	upLimit <- mean(offspring[x,]) + 2*sd(offspring[x,])
	if(any(offspring[x,]<downLimit)||any(offspring[x,]>upLimit)){
		result <- NULL
	}else{
		nrDistributions <- length(proportion)
		result <- rep(0,length(offspring[x,]))
		EM <- normalmixEM(sort(offspring[x,]), k=nrDistributions, maxrestarts=0, maxit = 100,fast=TRUE)
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
		if(checkMu.internal(offspring,EM,overlapInd)){
			result <- filterRow.internal(result, overlapInd, proportion, margin, genotypes)
		}else{
			result<- NULL
		}
	}
	invisible(result)
}

############################################################################################################
#									*** checkMu.internal ***
#
# DESCRIPTION:
#	checking if fitted normal distributions do not overlap
# 
# PARAMETERS:
# 	offspring - currently processed row
# 	EM - output of normalmixEM function
# 	overlapInd - how many individuals are allowed to be overlapping between distributions
# 
# OUTPUT:
#	boolean
#
############################################################################################################
checkMu.internal <- function(offspring,EM,overlapInd){
	for(i in 2:length(EM$mu)){
		up <- EM$mu[i]-EM$sigma[i]
		down <- EM$mu[i-1]+EM$sigma[i-1]
		if((up)<(down)){
			if(sum(offspring<down && offspring>up)){
				return(FALSE)
			}
		}		
	}
	return(TRUE)
}
