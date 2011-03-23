#################################################################################
#
# toGenotypes.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified March 2011
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
#           transformIndividual, zero, cEquals, cMore, cLess, checkExpression
#
#################################################################################

#toGenotypes: Function that chooses from the matrix only appropriate markers with specified rules

toGenotypes <- function(ril, treshold=0.01, overlapInd = 0, proportion = 50, margin = 5, verbose=FALSE, debugMode=0){
	#*******CHECKS*******
	if(proportion < 1 || proportion > 99) stop("Proportion is a percentage (1,99)")
	if(overlapInd < 0 || overlapInd > ncol(expressionMatrix)) stop("overlapInd is a number (0,lenght of the row).")
	if(margin < 0 || margin > proportion) stop("Margin is a percentage (0,proportion)")
	if(verbose && debugMode==1) cat("toGenotypes starting withour errors in checkpoint.\n")
	
	s <- proc.time()
	#*******SELECTING DIFFERENTIALLY EXPRESSED GENES*******
	s1 <- proc.time()
	ril <- selectDifferentiallyExpressed(ril,treshold,verbose,debugMode)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Selecting diff. expressed genes done in:",(e1-s1)[3],"seconds.\n")
	
	#*******CONVERTING CHILDREN PHENOTYPIC DATA TO GENOTYPES*******
	s1 <- proc.time()
	ril <- convertToGenotypes(ril, verbose, debugMode)
	ril <- filterGenotypes(ril, overlapInd, proportion, margin, verbose, debugMode)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Converting children phenotypic data to genotypes done in:",(e1-s2)[3],"seconds.\n")
	
	#*******SAVING CROSS OBJECT*******
	s1 <- proc.time()
	ril <- genotypesToCross(ril,verbose,debugMode)
	e1 <- proc.time()
	if(verbose && debugMode==2)cat("Converting children phenotypic data to genotypes done in:",(e1-s2)[3],"seconds.\n")
	
}

selectDifferentiallyExpressed( <- function(ril,treshold=0.01,verbose=FALSE,debugMode=0){
	s2<-proc.time()
	ril$parental$up$val <- ril$parental$phenotypes[which(ril$parental$RP$pval[1] < treshold),]
	ril$parental$up$means <- apply(ril$parental$up$val,1,mean)
	ril$parental$down$val <- ril$parental$phenotypes[which(ril$parental$RPpval[2] < treshold),]
	ril$parental$down$means <- apply(ril$parental$up$val,1,mean)
	ril$rils$up <- ril$rils$phenotypes[which(rownames(ril$rils$phenotypes) %in% rownames(ril$parental$up$val))]
	ril$rils$down <- ril$rils$phenotypes[which(rownames(ril$rils$phenotypes) %in% rownames(ril$parental$down$val))]
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Filtering data with treshold:",treshold,"done in:",(e2-s2)[3],"seconds.\n")
	invisible(ril)
}

convertToGenotypes <- function(ril,verbose=FALSE,debugMode=0){
	up <- lapply(rownames(ril$rils$up),splitRow,ril$rils$up,ril$parental$up)
	down <- lapply(rownames(ril$rils$down),splitRow,ril$rils$down,ril$parental$down)
	ril$rils$genotypes$simulated <- rbind(up,1-down)
	invisible(ril)
}

splitRow <- function(x,rils,parental){
	result <- rils[x,]
	splitVal <- parental$means[which(rownames(parental$val %in% x)),]
	result[which(rils[x,] > splitVal)] <- 0
	result[which(rils[x,] < splitVal)] <- 1
	result[which(rils[x,] == splitVal)] <- NA
	invisible(result)
}

filterGenotypes <- function(ril, overlapInd, proportion, margin, verbose=FALSE,debugMode=0){
	result <- apply(ril$rils$genotypes$simulated,1,filterRow,overlapInd, proportion, margin)
	ril$rils$genotypes$simulated <- ril$rils$genotypes$simulated[which(result==1),]
	invisible(ril)
}

filterRow <- function(genotypeRow, overlapInd, proportion, margin){
	if(sum(is.na(genotypeRow))!=0) return(0)
	above <- sum(genotypeRow==1)
	bellow <- sum(genotypeRow==0)
	if((above < (proportion+(margin/2))) && (above > (proportion-(margin/2)))){
		if((bellow < ((100-proportion)+(margin/2))) && (bellow > ((100-proportion)-(margin/2)))){
			return(1)
		}
	}
	return(0)
}
