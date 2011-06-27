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
#           convertToGenotypes.internal, splitRow.internal, filterGenotypes.internal, filterRow.internal
#			sortMap.internal, majorityRule.internal, mergeChromosomes.internal, splitRowSubEM.internal
#			checkMu.internal
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
# 	... - Parameters passed to formLinkageGroups.
# 	verbose - Be verbose
# 	debugMode - 1: Print our checks, 2: print additional time information
# 
# OUTPUT:
#	object of class cross
#
############################################################################################################
toGenotypes <- function(population, genotype=c("simulated","real"), orderUsing=c("map_genetic","map_physical"), splitMethod=c("EM","mean"),treshold=0.01, overlapInd = 0, proportion = list(c(50,50),15), margin = 15, removalFUN = NULL, verbose=FALSE, debugMode=0,...){
	#*******CHECKS*******
	if(missing(orderUsing)) orderUsing <- NULL
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
	inListCheck.internal(orderUsing,"orderUsing",c("map_genetic","map_physical"))
	
	
	#*******CONVERTING CHILDREN PHENOTYPIC DATA TO GENOTYPES*******
	if(genotype=="simulated"){
		s1 <- proc.time()
		population <- convertToGenotypes.internal(population, splitMethod, treshold, overlapInd, proportion, margin, verbose, debugMode)
		e1 <- proc.time()
		if(verbose && debugMode==2)cat("Converting phenotypes to genotypes done in:",(e1-s1)[3],"seconds.\n")
	}
	
	
	#*******SAVING CROSS OBJECT*******
	s1 <- proc.time()
	cross <- genotypesToCross.internal(population,genotype=genotype,orderUsing=orderUsing,verbose=verbose,debugMode=debugMode)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Creating cross object done in:",(e1-s1)[3],"seconds.\n")
	
	#*******ENHANCING CROSS OBJECT*******
	#if(use!="map_genetic"&&use!="map_physical"){
		### FormLinkage groups
	#	cross <- invisible(formLinkageGroups(cross,reorgMarkers=TRUE,verbose=verbose,...))
		
		### remove shitty chromosomes
	#	if(!is.null(numberOfChromosomes))	cross <- removeChromosomes.internal(cross,numberOfChromosomes)
		### saving as separated object, beacause orderMarkers will remove it from cross object
	#	removed <- cross$rmv
		
		### Order markers
		#cross <- orderMarkers(cross, use.ripple=TRUE, verbose=verbose)
		
		### Adding real maps		
	#	if(!(is.null(population$maps$physical))){
	#		population <- sortMap.internal(population)
	#		cross$maps$physical <- population$maps$physical
			### Majority rule used to order linkage groups
	#		cross <- segregateChromosomes.internal(cross)
	#	}
		
		### adding info about removed markers
	#	cross$rmv <- removed
	#}
	
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
convertToGenotypes.internal <- function(population, splitMethod, treshold, overlapInd, proportion, margin, verbose=FALSE, debugMode=0){
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
	upRils <- population$offspring$phenotypes[rownames(upParental),]
	### down-regulated
	downNotNull <- which(population$founders$RP$pval[2] > 0)
	downBelowTreshold <- which(population$founders$RP$pval[2] < treshold)
	downSelected <- downBelowTreshold[which(downBelowTreshold%in%downNotNull)]
	downParental <- population$founders$phenotypes[downSelected,]
	downRils <- population$offspring$phenotypes[rownames(downParental),]
	
	### checking if anything is selected and if yes - processing
	if(!(is.null(dim(upRils)))){
		if(!(is.null(dim(downRils)))){
			# best situation
			if(verbose) cat("Selected ",nrow(upRils),"upregulated markers and ",nrow(downRils),"downregulated markers.\n")
			cur <- splitRow.internal(downRils, downParental, splitMethod, overlapInd, proportion, margin, population$founders$groups, 0)
			output <- rbind(output,cur[[1]])
			markerNames <- c(markerNames,cur[[2]])
		}else{
			if(verbose) cat("Selected ",nrow(upRils),"upregulated markers.\n")
		}
		cur <- splitRow.internal(upRils, upParental, splitMethod, overlapInd, proportion, margin, population$founders$groups, 1)
		output <- rbind(output,cur[[1]])
		markerNames <- c(markerNames,cur[[2]])
	}else{
		if(!(is.null(dim(downRils)))){
			if(verbose) cat("Selected ",nrow(downRils),"downregulated markers.\n")
			cur <- splitRow.internal(downRils, downParental, splitMethod, overlapInd, proportion, margin, population$founders$groups, 0)
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
#									*** splitRow.internal ***
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
############################################################################################################
splitRow.internal <- function(offspring, founders, splitMethod, overlapInd, proportion, margin, groupLabels, up){
	output <- NULL
	markerNames <- NULL
	cat("population:",nrow(offspring),"founders:",nrow(founders),"\n")
	for(x in rownames(offspring)){
		if(splitMethod=="mean"){
			cur <- splitRowSub.internal(x, offspring, founders, overlapInd, proportion, margin, groupLabels, up)
		}else if(splitMethod=="EM"){
			cur <- splitRowSubEM.internal(x, offspring, founders, overlapInd, proportion, margin, groupLabels, up)
		}
		if(!(is.null(cur))){
			output <- rbind(output,cur)
			markerNames <- c(markerNames,x)
		}
	}
	invisible(list(output,markerNames))
}

############################################################################################################
#splitRowSub.internal: subfunction of splitRow.internal, splitting one row
# 
# x - name of currently processed row
# offspring - matrix of up/down regulated genes in offspring
# founders - matrix of up/down regulated genes in parents
# overlapInd - Number of individuals that are allowed in the overlap
# proportion - Proportion of individuals expected to carrying a certain genotype 
# margin - Proportion is allowed to varry between this margin (2 sided)
# groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
# up - 1 - genes up 0 - down regulated
#
############################################################################################################
splitRowSub.internal <- function(x, offspring, founders, overlapInd, proportion, margin, groupLabels, up=1){
	### initialization
	result <- rep(0,length(offspring[x,]))
	
	### splitting
	if(length(proportion)==2){
		### splitting into 0/1 genotypes
		if(up==1){
			genotypes <- c(0,1)
		}else if(up==0){
			genotypes <- c(1,0)
		}
		A <- founders[x,which(groupLabels==0)]
		B <- founders[x,which(groupLabels==1)]
		splitVal <- mean(mean(A,na.rm=TRUE),mean(B,na.rm=TRUE))
		offspring[x,which(is.na(offspring[x,]))] <- splitVal
		a <- offspring[x,] > splitVal
		b <- offspring[x,] < splitVal
		result[which(a)] <- genotypes[1]
		result[which(b)] <- genotypes[2]
		result[which(offspring[x,] == splitVal)] <- NA
		result <- filterRow.internal(result, overlapInd, proportion, margin, genotypes)
	}else if(length(proportion)==3){
		### splitting into 0/1/2 genotypes
		if(up==1){
			genotypes <- c(0,1,2)
		}else if(up==0){
			genotypes <- c(2,1,0)
		}
		meanVal <- mean(offspring[x,],na.rm=TRUE)
		sdVal <- sd(offspring[x,],na.rm=TRUE)
		subSplitValDown <- meanVal - sdVal
		subSplitValUp <- meanVal + sdVal
		result[which(offspring[x,] < subSplitValDown )] <- genotypes[1]
		result[which(offspring[x,] > subSplitValUp )] <- genotypes[3]
		result[which((offspring[x,] > subSplitValDown )&&(offspring[x,] < subSplitValUp ))] <- genotypes[2]
		result[which(offspring[x,] == meanVal)] <- NA
		result <- filterRow.internal(result, overlapInd, proportion, margin, genotypes)
	}
	
	invisible(result)
}

############################################################################################################
#filterRow.internal : removing from genotypic matrix genes that are no passing specified requirments
# 
# population - Ril type object, must contain founders phenotypic data.
# overlapInd - Number of individuals that are allowed in the overlap
# proportion - Proportion of individuals expected to carrying a certain genotype 
# margin - Proportion is allowed to varry between this margin (2 sided)
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information

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
	if(filterRowSub.internal(result,overlapInd,proportion, margin, genotypes)){
		invisible(result)
	}else if(filterRowSub.internal(result,overlapInd,proportion, margin, genotypes)){
		invisible(result)
	}else{
		return(NULL)
	}
}

############################################################################################################
#filterRowSub.internal : subfunction of filterGenotypes.internal, filtering one row
# 
# genotypeRow - currently processed row
# overlapInd - Number of individuals that are allowed in the overlap
# proportion - Proportion of individuals expected to carrying a certain genotype 
# margin - Proportion is allowed to varry between this margin (2 sided)
#
############################################################################################################
filterRowSub.internal <- function(result, overlapInd, proportion, margin, genotypes){
	for(i in 1:length(proportion)){
		cur <- sum(result==genotypes[i])/length(result) * 100
		if(is.na(cur)) return(FALSE)
		upLimit <- proportion[i] + margin/2
		downLimit <- proportion[i] - margin/2
		#cat(result,"\n",cur,upLimit,downLimit,"\n")
		if(!((cur < upLimit)&&(cur > downLimit))){
			return(FALSE)
		}
	}
	return(TRUE)
}

############################################################################################################
#sortMap.internal: subfunction of filterGenotypes.internal, filtering one row
# 
# population - Ril type object
#
############################################################################################################
sortMap.internal <- function(population){
	genes <- population$offspring$map
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
	population$offspring$map <- list(result,lengths[-length(lengths)])
	invisible(population)
}

############################################################################################################
#segragateChromosomes.internal  - ordering chromosomes using physical map
# 
# cross - object of R/qtl cross type
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
#majorityRule.internal - subfunction of segragateChromosomes.internal, returns matrix showing for every
# reco map chromosome from which physicall map chromosome majority of markers comes
# 
# cross - object of R/qtl cross type
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
#mergeChromosomes.internal - subfunction of segragateChromosomes.internal, merging multiple chromosomes into
# one
# 
# cross - object of R/qtl cross type
# chromosomes - chromosomes to be merged
# name - chromosome
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

############################################################################################################
#splitRowSubEM.internal - subfunction of segragateChromosomes.internal, merging multiple chromosomes into
# one
# 
# cross - object of R/qtl cross type
# chromosomes - chromosomes to be merged
# name - chromosome
#
############################################################################################################
splitRowSubEM.internal <- function(x, offspring, founders, overlapInd, proportion, margin, groupLabels, up=1){
	### initialization
	print(x)
	nrDistributions <- length(proportion)
	result <- rep(0,length(offspring[x,]))
	EM <- cat(normalmixEM(sort(offspring[x,]), k=nrDistributions, maxrestarts=0, maxit = 100,fast=TRUE),file=NULL)
	if(up==1){
		genotypes <- c(0:(nrDistributions-1))
	}else if(up==0){
		genotypes <- c((nrDistributions-1):0)
	}
	len <- vector(mode="numeric",length=nrDistributions)
	for(i in 1:nrDistributions){
		len[i]<-length(offspring[x,])*EM$lambda[i]
		startVal <- sum(len[1:i-1])
		offspring[x,which(offspring[x,] %in% sort(offspring[x,])[startVal:(startVal+len[i])])] <- genotypes[i]
	}
	#result <- checkMu.internal(EM)
	if(checkMu.internal(EM)){
		result <- filterRow.internal(result, overlapInd, proportion, margin, genotypes)
	}else{
		result<- NULL
	}
	invisible(result)
}

############################################################################################################
#checkMu.internal: checking if fitted normal distributions do not overlap
# 
# EM - output of normalmixEM function
#
############################################################################################################
checkMu.internal <- function(EM){
	for(i in 2:length(EM$mu)){
		if((EM$mu[i]-EM$sigma[i])<(EM$mu[i-1]+EM$sigma[i-1])) return(FALSE)
	}
	return(TRUE)
}
