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

toGenotypes <- function(ril, treshold=0.01,verbose=FALSE,debugMode=0){
	ril <- selectDifferentiallyExpressed(ril)
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