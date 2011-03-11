#####################################################################
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
#####################################################################

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
toGenotypes <- function(expressionMatrix, splitFUN = zero, overlapInd = 0, proportion = 50, margin = 5, genotypes = c(0,1), verbose=FALSE, debugMode=0){
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
	
	if(verbose) cat("Done, appriopiateMarkers took:",(e-s)[3],"seconds.\n")
	
	invisible(genotypeMatrix)
}

#transformIndividual: 
# x
# r
# genotypes
transformIndividual <- function(x,r,genotypes){
	results <- x
	results[which(x  > r)] <- genotypes[1]
	results[which(x  < r)] <- genotypes[2]
	results[which(x  == r)] <- NA
	results
}

#zero: 
# x
zero <- function(x){
	return(0)
}

#cEquals
# counts how many times x eqauls splitFUN
cEquals <- function(x,splitFUN){
	sum(x==splitFUN(x))
}

#cLess
# counts how many times x is less than splitFUN
cLess <- function(x,splitFUN){
	sum(x>splitFUN(x))
}

#counts how many times x is more than splitFUN
cMore <- function(x,splitFUN){
	sum(x<splitFUN(x))
}

#check - checks if x meets specified requirments:
#splitFUN -> function used to split values
#overlapInd - how many individuals could be overalpping -> how many elememnts of x are eqaul to splitFUN
#margin_range - we assume that nr of less and more should be the same, at least difference should be less than margin_range
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

#test function, doesn't need documenation
test.toGenotypes <- function(){
	data(expressionData)
	expressionMatrix <- as.matrix(read.table("Expression_BrassicaRapa_10chr2.txt",sep=""))
	genotypes <- toGenotypes(expressionMatrix, margin=0.5,genotypes=c(1,0),overlapInd=0, verbose=TRUE, debugMode=2)
	
	#Checks
	if(sum(which(genotypes[1,1]!=1)))	stop("Element not equal to expected\n")
	if(sum(which(genotypes[100,10]!=1)))	stop("Element not equal to expected\n")
	if(sum(which(genotypes[1026,20]!=1)))	stop("Element not equal to expected\n")
	if(sum(dim(genotypes))!=2063)	stop("Wrong dimensions\n")
	
	invisible(genotypes)
}
