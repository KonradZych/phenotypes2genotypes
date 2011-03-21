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
#
# expressionMatrix -> columns: individuals, rows: markers
# proportion -> Proportion of individuals expected to carrying a certain genotype
# margin -> Proportion is allowed to varry between this margin (2 sided)
# splitFUN -> function used to split values
# overlapInd -> number of individuals that are allowed to overlap (have splitFUN)
# genotypes -> User defined genotypes for the output matrix
# verbose standard
# debugmode standard ->1 Print our checks, 2 print additional time information
toGenotypes <- function(expressionMatrix,  splitFUN = zero, overlapInd = 0, proportion = 50, margin = 5, genotypes = c(0,1), verbose=FALSE, debugMode=0){
	s <- proc.time()

	if(proportion < 1 || proportion > 99) stop("Proportion is a percentage (1,99)")
	if(overlapInd < 0 || overlapInd > ncol(expressionMatrix)) stop("overlapInd is a number (0,lenght of the row).")
	if(margin < 0 || margin > proportion) stop("Margin is a percentage (0,proportion)")
	if(verbose && debugMode==1) cat("toGenotypes starting withour errors in checkpoint.\n")

	#Selection of the probes matching to the specified parameters
	suitedRows <- apply(expressionMatrix,1,checkExpression,splitFUN,overlapInd, proportion, margin)
	expressionMatrix <- expressionMatrix[which(suitedRows),]

	ep <- proc.time()
	if(verbose && debugMode==2) cat("Selected proper probes, took:",(ep-s)[3],"seconds. Creating genotype matrix.\n")

	#Transform numeric values to genotypes
	r <- apply(expressionMatrix,1,splitFUN)
	genotypeMatrix <- apply(expressionMatrix,2,transformIndividual,r,genotypes)

	eg <- proc.time()
	if(verbose && debugMode==2) cat("Created genotype matrix, took:",(eg-ep)[3],"seconds.\n")

	e<-proc.time()
	if(verbose) cat("toGenotypes finished in",(e-s)[3],"seconds.\n")
	
	invisible(genotypeMatrix)
}

#transformIndividual: spliting column form expressionMatrix using vector of spliting values
# x - column of expressionMatrix, containing data 
# r - vector containing split values (result of applying specified function throughtoutthe rows of expressionMatrix)
# genotypes - list of 2 values, specifing how genotypes should be named
transformIndividual <- function(x,r,genotypes){
	results <- x
	results[which(x  > r)] <- genotypes[1]
	results[which(x  < r)] <- genotypes[2]
	results[which(x  == r)] <- NA
	results
}

parentalSplit <- function(x,expressionChildren,expressionParental,verbose=FALSE, debugMode=0){
	if(verbose && debugMode==1)if(x%%100==0)cat("parentalSplit starting withour errors in checkpoint for row:",x,"\n")
	s <- proc.time()
	expressionParentalRow <- expressionParental[which(rownames(expressionParental) %in% rownames(expressionChildren)[x]),]
	genotypeMatrixRow <- lapply(expressionChildren[x,],parentalSplitSub,expressionParentalRow)
	e <- proc.time()
	if(verbose && debugMode==2)if(x%%100==0)cat("parentalSplit for row:",x,"done in:",(e-s)[3],"seconds.\n")
	invisible(genotypeMatrixRow)
}

parentalSplitSub <- function(expressionChildrenElement,expressionParentalRow){
	distance1 <- abs(expressionChildrenElement-expressionParentalRow[1])
	distance2 <- abs(expressionChildrenElement-expressionParentalRow[2])
	if(distance1<=distance2){
		genotypeMatrixElement <- 0
	}else{
		genotypeMatrixElement <- 1
	}
	invisible(genotypeMatrixElement)
}

#zero: function returning 0, to use with toGenotypes, when splitting value is exactly 0
# x - input of any type, when used with toGenotypes - row of an expressionMatrix
zero <- function(x){
	return(0)
}

#cEquals
# counts how many elements of x eqaul result of splitFUN(x)
cEquals <- function(x,splitFUN){
	sum(x==splitFUN(x))
}

#cLess
# counts how many elements of x are less than result of splitFUN(x)
cLess <- function(x,splitFUN){
	sum(x>splitFUN(x))
}

#cMore
# counts how many elements of x are more than result of splitFUN(x)
cMore <- function(x,splitFUN){
	sum(x<splitFUN(x))
}

#check - checks if x meets specified requirments:
# splitFUN -> function used to split values
# overlapInd - how many individuals could be overalpping -> how many elememnts of x are eqaul to result of splitFUN(x)
# margin_range - we assume that nr of less and more should be the same, at least difference should be less than margin_range
checkExpression <- function(x, splitFUN, overlapInd, proportion, margin){
	r <- FALSE
	if(cEquals(x,splitFUN) <= overlapInd){
		above <- cMore(x,splitFUN) / length(x) * 100
		bellow <- cLess(x,splitFUN) / length(x) * 100
		if((above < (proportion+(margin/2))) && (above > (proportion-(margin/2)))){
		if((bellow < ((100-proportion)+(margin/2))) && (bellow > ((100-proportion)-(margin/2)))){
			r <- TRUE
		}
		}
	}
	r
}
