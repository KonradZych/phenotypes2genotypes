############################################################################################################
#
# toGenotypes.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified April 2011
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
#			sortMap.internal 
#
############################################################################################################

############################################################################################################
#toGenotypes: Function that chooses from the matrix only appropriate markers with specified rules
# 
# ril - Ril type object, must contain parental phenotypic data.
# use - Which genotypic matrix should be saved to file, real - supported by user and read from file, 
#	simulated - made by toGenotypes, ap - simulated data orderd using gff map
# treshold - If Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.
# overlapInd - Number of individuals that are allowed in the overlap
# proportion - Proportion of individuals expected to carrying a certain genotype 
# margin - Proportion is allowed to varry between this margin (2 sided)
# minChrLength -if maximal distance between the markers in the chromosome is lower than this value,
#	whole chromosome will be dropped
# ... - Parameters passed to formLinkageGroups.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
toGenotypes <- function(ril, use=c("real","simulated","map"), treshold=0.01, overlapInd = 0, proportion = 50, margin = 15, minChrLength = 0, verbose=FALSE, debugMode=0,...){
	#*******CHECKS*******
	require(qtl)
	if(proportion < 1 || proportion > 99) stop("Proportion is a percentage (1,99)")
	#if(overlapInd < 0 || overlapInd > ncol(expressionMatrix)) stop("overlapInd is a number (0,lenght of the row).")
	if(margin < 0 || margin > proportion) stop("Margin is a percentage (0,proportion)")
	if(verbose && debugMode==1) cat("toGenotypes starting withour errors in checkpoint.\n")
	
	#*******CONVERTING CHILDREN PHENOTYPIC DATA TO GENOTYPES*******
	s1 <- proc.time()
	ril <- convertToGenotypes.internal(ril, treshold, verbose, debugMode)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Converting phenotypes to genotypes done in:",(e1-s1)[3],"seconds.\n")
	
	#*******FILTERING GENOTYPES*******
	s1 <- proc.time()
	ril <- filterGenotypes.internal(ril, overlapInd, proportion, margin, verbose, debugMode)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Selecting markers using specified parameters done in:",(e1-s1)[3],"seconds.\n")
	
	#*******SAVING CROSS OBJECT*******
	s1 <- proc.time()
	cross <- genotypesToCross.internal(ril,use=use,verbose=verbose,debugMode=debugMode)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Creating cross object done in:",(e1-s1)[3],"seconds.\n")
	
	#*******ENHANCING CROSS OBJECT*******
	if(use!="map"){
		#FormLinkage groups
		cross <- invisible(formLinkageGroups(cross,reorgMarkers=TRUE,verbose=verbose,...))
		### remove shitty chromosomes
		cross <- removeChromosomes.internal(cross,minChrLength)
		removed <- cross$rmv
		#Order markers - ripple in our case often produces errors, so turning it off by default
		cross <- orderMarkers(cross, use.ripple=TRUE, verbose=verbose)
		#Adding real maps
		
		if(!(is.null(ril$rils$map))){
			ril <- sortMap.internal(ril)
			cross$maps$physical <- ril$rils$map
		}
		cross$rmv <- removed
		#Majority rule
		###********IMPORTANT BUT NOT YET OPERATIVE!!!********
		#cross <- segragateChromosomes.internal(cross)
	}
	
	#*******RETURNING CROSS OBJECT*******
	invisible(cross)
}

############################################################################################################
#convertToGenotypes.internal: function splitting differentially expressed markers into two genotypes
# 
# ril - Ril type object, must contain parental phenotypic data.
# treshold - If Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information 
#
############################################################################################################
convertToGenotypes.internal <- function(ril,treshold,verbose=FALSE,debugMode=0){
	if(verbose && debugMode==1) cat("convertToGenotypes starting.\n")
	m <- NULL
	upParental <- ril$parental$phenotypes[which(ril$parental$RP$pval[1] < treshold),]
	downParental <- ril$parental$phenotypes[which(ril$parental$RP$pval[2] < treshold),]
	upRils <- ril$rils$phenotypes[which(rownames(ril$rils$phenotypes) %in% rownames(upParental)),]
	downRils <- ril$rils$phenotypes[which(rownames(ril$rils$phenotypes) %in% rownames(downParental)),]
	for(x in rownames(upRils)){
		m <- rbind(m,splitRow.internal(x,upRils,upParental,ril$parental$groups,c(0,1)))
	}
	c<-0
	for(x in rownames(downRils)){
		m <- rbind(m,splitRow.internal(x,downRils,downParental,ril$parental$groups,c(1,0)))
	}
	ril$rils$genotypes$simulated <- m
	colnames(ril$rils$genotypes$simulated) <- colnames(upRils)
	rownames(ril$rils$genotypes$simulated) <- c(rownames(upRils),rownames(downRils))
	invisible(ril)
}

############################################################################################################
#splitRow.internal: subfunction of convertToGenotypes.internal, splitting one row
# 
# x - name of currently processed row
# rils - matrix of up/down regulated genes in rils
# parental - matrix of up/down regulated genes in parents
# genotypes - values genotypic matrix will be filled with
#
############################################################################################################
splitRow.internal <- function(x,rils,parental,groupLabels,genotypes){
	result <- rep(0,length(rils[x,]))
	A <- parental[which(rownames(parental) == x),which(groupLabels==0)]
	B <- parental[which(rownames(parental) == x),which(groupLabels==1)]
	splitVal <- mean(mean(A,na.rm=TRUE),mean(B,na.rm=TRUE))
	result[which(rils[x,] > splitVal)] <- genotypes[1]
	result[which(rils[x,] < splitVal)] <- genotypes[2]
	result[which(rils[x,] == splitVal)] <- NA
	invisible(result)
}

############################################################################################################
#filterGenotypes.internal: removing from genotypic matrix genes that are no passing specified requirments
# 
# ril - Ril type object, must contain parental phenotypic data.
# overlapInd - Number of individuals that are allowed in the overlap
# proportion - Proportion of individuals expected to carrying a certain genotype 
# margin - Proportion is allowed to varry between this margin (2 sided)
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information

############################################################################################################
filterGenotypes.internal <- function(ril, overlapInd=0, proportion=50, margin=5, verbose=FALSE,debugMode=0){
	if(verbose && debugMode==1) cat("filterGenotypes starting.\n")
	result <- apply(ril$rils$genotypes$simulated,1,filterRow.internal,overlapInd=overlapInd,proportion=proportion, margin=margin)
	ril$rils$genotypes$simulated <- ril$rils$genotypes$simulated[which(result==1),]
	invisible(ril)
}

############################################################################################################
#filterRow.internal: subfunction of filterGenotypes.internal, filtering one row
# 
# genotypeRow - currently processed row
# overlapInd - Number of individuals that are allowed in the overlap
# proportion - Proportion of individuals expected to carrying a certain genotype 
# margin - Proportion is allowed to varry between this margin (2 sided)
#
############################################################################################################
filterRow.internal <- function(genotypeRow, overlapInd, proportion, margin){
	if(sum(is.na(genotypeRow))>overlapInd) return(0)
	above <- sum(genotypeRow==1)/length(genotypeRow) * 100
	below <- sum(genotypeRow==0)/length(genotypeRow) * 100
	if((above < (proportion+margin/2)) && (above > (proportion-margin/2))){
		if((below > 100-(proportion+margin/2)) && (below < 100-(proportion-margin/2))){
			return(1)
		}
	} else if((below < (proportion+margin/2)) && (below > (proportion-margin/2))){
		if((above > 100-(proportion+margin/2)) && (above < 100-(proportion-margin/2))){
			return(1)
		}
	}
	return(0)
}

############################################################################################################
#sortMap.internal: subfunction of filterGenotypes.internal, filtering one row
# 
# ril - Ril type object
#
############################################################################################################
sortMap.internal <- function(ril){
	genes <- ril$rils$map
	result <- NULL
	nchr <- length(table(genes[,1]))
	lengths <- vector(mode="numeric",length=nchr+1)
	lengths[1] <- 0
	for(i in 1:nchr){
		current <- genes[which(genes[,1]==i),]
		lengths[i+1] <- max(current[,2]) + lengths[i] + 20000
		current <- current[names(sort(current[,2])),]
		current[,2] <- current[,2] + lengths[i]
		result <- rbind(result,current)
	}
	ril$rils$map <- list(result,lengths[-length(lengths)])
	invisible(ril)
}
