#################################################################################
#
# parentalRoutine.R
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
# Contains: parentalRoutine 
# 				readParentalExpression, rankParentalExpression, filterParentalExpression
#
#################################################################################

readFiles <- function(rils="children",parental="parental",verbose=FALSE,debugMode=0){
	s <- proc.time()
	ril <- NULL
	
	filename <- paste(rils,"_phenotypes.txt",sep="")
	if(file.exists(filename)){
		ril$rils$phenotypes <- readFile(filename,verbose,debugMode)
		if(verbose) cat("Found phenotypic file for rils:",filename,"and stored it in ril$rils$phenotypes\n")
	}else{
		stop("There is no phenotypic file for rils: ",filename," this file is essentiall, you have to provide it\n")
	}
	
	filename <- paste(rils,"_genotypes.txt",sep="")
	if(file.exists(filename)){
		ril$rils$genotypes$read <- readFile(filename,verbose,debugMode)
		if(verbose) cat("Found genotypic file for rils:",filename,"and stored it in ril$rils$genotypes\n")
	}else{
		if(verbose) cat("There is no genotypic file for rils:",filename,"genotypic data for rils will be simulated\n")
	}
	
	filename <- paste(parental,"_phenotypes.txt",sep="")
	if(file.exists(filename)){
		ril$parental$phenotypes <- readFile(filename,verbose,debugMode)
		if(verbose) cat("Found phenotypic file for parents:",filename,"and stored it in ril$parental$phenotypes\n")
	}else{
		if(verbose) cat("There is no phenotypic file for parents:",filename,"further processing will take place without taking into account parental data\n")
	}
	
	e <- proc.time()
	return(ril)
}

readFile <- function(filename,verbose=FALSE,debugMode=0){
	#CRUCIAL CHECKS
	s1<-proc.time()
	expressionMatrix <- as.matrix(read.table(filename,sep=""))
	e1<-proc.time()
	if(verbose && debugMode==2)cat("Reading expression file:",filename,"done in:",(e1-s1)[3],"seconds.\n")
	invisible(expressionMatrix)
}

