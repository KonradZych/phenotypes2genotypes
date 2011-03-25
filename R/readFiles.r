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
# Contains: readFiles 
# 				readFile.internal, mapMarkers.internal
#
#################################################################################

readFiles <- function(rils="children",parental="parental",sep="",verbose=FALSE,debugMode=0){
	#**********INITIALIZING FUNCTION*************
	s <- proc.time()
	if(verbose && debugMode==1) cat("readFiles starting.\n")
	invisible(require(iqtl))
	if(verbose) cat("Loaded required libraries.\n")
	ril <- NULL
	
	#**********READING CHILDREN PHENOTYPIC DATA*************
	filename <- paste(rils,"_phenotypes.txt",sep="")
	if(file.exists(filename)){
		if(verbose) cat("Found phenotypic file for rils:",filename,"and will store  it in ril$rils$phenotypes\n")
		ril$rils$phenotypes <- readFile.internal(filename,sep,verbose,debugMode)
	}else{
		stop("There is no phenotypic file for rils: ",filename," this file is essentiall, you have to provide it\n")
	}
	
	#**********READING CHILDREN GENOTYPIC DATA*************
	filename <- paste(rils,"_genotypes.txt",sep="")
	if(file.exists(filename)){
		if(verbose) cat("Found genotypic file for rils:",filename,"and will store  it in ril$rils$genotypes$read\n")
		ril$rils$genotypes$read <- readFile.internal(filename,sep,verbose,debugMode)
	}else{
		cat("WARNING: There is no genotypic file for rils:",filename,"genotypic data for rils will be simulated\n")
	}
	
	#**********READING PARENTAL PHENOTYPIC DATA*************
	filename <- paste(parental,"_phenotypes.txt",sep="")
	if(file.exists(filename)){
		if(verbose) cat("Found phenotypic file for parents:",filename,"and will store it in ril$parental$phenotypes\n")
		ril$parental$phenotypes <- readFile.internal(filename,sep,verbose,debugMode)
		#removing from parental probes that are not in children:
		ril$parental$phenotypes <- mapMarkers.internal(ril$parental$phenotypes,ril$rils$phenotypes ,mapMode=1,verbose=verbose)
	}else{
		cat("WARNING: There is no phenotypic file for parents:",filename,"further processing will take place without taking into account parental data\n")
	}
	
	#**********FINALIZING FUNCTION*************
	e <- proc.time()
	if(verbose) cat("readFiles done in",(e-s)[3],"seconds.\n")
	invisible(ril)
}

readFile.internal <- function(filename,sep="",verbose=FALSE,debugMode=0){
	if(!file.exists(filename)) stop("File: ",filename,"doesn't exist.\n")
	s1<-proc.time()
	currentFile <- read.table(filename,sep=sep)
	currentFile <- as.matrix(currentFile)
	if(!is.numeric(currentFile)) stop("Not all elements of the file:",filename," are numeric, check help files for rigth input file format. \n")
	if(verbose) cat("File:",filename,"was loaded into R enviroment containing",nrow(currentFile),"markers and",ncol(currentFile),"individuals.\n")
	e1<-proc.time()
	if(verbose && debugMode==2)cat("Reading file:",filename,"done in:",(e1-s1)[3],"seconds.\n")
	invisible(currentFile)
}

mapMarkers.internal <- function(expressionMatrix1, expressionMatrix2, mapMode=2,verbose=FALSE,debugMode=0){
	if(mapMode==1) {
		nrRows <- nrow(expressionMatrix1)
		#warnings when names are mismatching
		if(verbose && debugMode==2)if(nrRows!=nrow(expressionMatrix2)){
			cat("Following markers will be removed:\n")
			cat(paste(rownames(expressionMatrix1)[which(!(rownames(expressionMatrix1) %in% rownames(expressionMatrix2)))],"\n"))
		}
		#mapping itself
		expressionMatrix1 <- expressionMatrix1[which(rownames(expressionMatrix1) %in% rownames(expressionMatrix2)),]
		if(verbose) cat("Because of names mismatch,",nrRows-nrow(expressionMatrix1),"markers were removed, run function with verbose=T debugMode=2 to print their names out.\n")
	}
	else if(mapMode==2){
		nrCols <- ncol(expressionMatrix1)
		#warnings when names are mismatching
		if(verbose && debugMode==2)if(nrCols!=ncol(expressionMatrix2)){
			cat("Following individuals will be removed:\n")
			paste(colnames(expressionMatrix1)[which(!(colnames(expressionMatrix1) %in% colnames(expressionMatrix2)))],"\n")
		}
		#mapping itself
		expressionMatrix1 <- expressionMatrix1[,which(colnames(expressionMatrix1) %in% colnames(expressionMatrix2))]
		if(verbose) cat("Because of names mismatch,",nrCols-ncol(expressionMatrix1),"individuals were removed, run function with verbose=T debugMode=2 to print their names out.\n")
	}
	invisible(expressionMatrix1)
}

