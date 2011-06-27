#############################################################################################
#
# readFiles.R
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
# Contains: readFiles 
# 				mapMarkers.internal, gffParser, probesLocation.internal
#				correctRowGff.internal, correctRowLoc.internal
#
#############################################################################################

############################################################################################################
#									*** readFiles ***
#
# DESCRIPTION:
#	reads geno/phenotypic files into R environment into special object.
# 
# PARAMETERS:
# 	offspring - Core used to specify names of children phenotypic ("offspring_phenotypes.txt") and genotypic ("offspring_genotypes.txt") files.
# 	founders - Core used to specify names of founders phenotypic ("founders_phenotypes.txt") file.
# 	map - Core used to specify names of genetic ("map_genetic.txt") and physical ("map_physical.txt") map files.
# 	sep - Separator of values in files. Passed directly to read.table, so "" is a wildcard meaning whitespace.
# 	verbose - Be verbose
# 	debugMode - 1: Print our checks, 2: print additional time information
#
# OUTPUT:
#	object of class population 
#
############################################################################################################
readFiles <- function(offspring="offspring",founders="founders",map="maps",sep="\t",verbose=FALSE,debugMode=0){
	#**********INITIALIZING FUNCTION*************
	s <- proc.time()
	if(verbose && debugMode==1) cat("readFiles starting.\n")
	population <- NULL
	
	#**********READING CHILDREN PHENOTYPIC DATA*************
	filename <- paste(offspring,"_phenotypes.txt",sep="")
	if(file.exists(filename)){
		if(verbose) cat("Found phenotypic file for offspring:",filename,"and will store  it in population$offspring$phenotypes\n")
		offspring_phenotypes <- read.table(filename,sep=sep,header=TRUE)
		offspring_phenotypes <- as.matrix(offspring_phenotypes)
		population <- createPopulation(offspring_phenotypes,no.warn=TRUE)
		doCleanUp.internal()
	}else{
		stop("No phenotype file for offspring: ",filename," this file is essentiall, you have to provide it\n")
	}
	
	#**********READING CHILDREN GENOTYPIC DATA*************
	filename <- paste(offspring,"_genotypes.txt",sep="")
	if(file.exists(filename)){
		if(verbose) cat("Found genotypic file for offspring:",filename,"and will store  it in population$offspring$genotypes$read\n")
		offspring_genotypess <- read.table(filename,sep=sep,header=TRUE)
		offspring_genotypess <- as.matrix(offspring_genotypess)
		population <- intoPopulation(population, offspring_genotypess, "offspring$genotypes")
		doCleanUp.internal()
	}else{
		if(verbose)cat("No genotypic file for offspring:",filename,"genotypic data for offspring will be simulated\n")
	}
	
	#**********READING PARENTAL PHENOTYPIC DATA*************
	filename <- paste(founders,"_phenotype.txt",sep="")
	if(file.exists(filename)){
		if(verbose) cat("Found phenotypic file for parents:",filename,"and will store it in population$founders$phenotypes\n")
		founders <- read.table(filename,sep=sep,header=TRUE)
		founders <- as.matrix(founders)
		population <- intoPopulation(population, founders, "founders")
		#removing from founders probes that are not in children:
		population$founders$phenotypes <- mapMarkers.internal(population$founders$phenotypes,population$offspring$phenotypes, mapMode=1, verbose=verbose)
		doCleanUp.internal()
	}else{
		warning("No phenotype file for parents: ",filename,". Strongly recommend to supply this data.\n")
	}
	
	#**********READING GENETIC MAP*************
	filename <- paste(map,"_genetic.txt",sep="")
	if(file.exists(filename)){
		if(verbose) cat("Found genotypic file for offspring:",filename,"and will store  it in population$offspring$genotypes$read\n")
		maps_genetic <- read.table(filename,sep=sep,row.names=1,header=FALSE)
		maps_genetic <- as.matrix(maps_genetic)
		population <- intoPopulation(population, maps_genetic, "maps$genetic")
		doCleanUp.internal()
	}else{
		if(verbose)cat("No genetic map file:",filename,".\n")
	}
	
	#**********READING PHYSICAL MAP*************
	filename <- paste(map,"_physical.txt",sep="")
	if(file.exists(filename)){
		if(verbose) cat("Found genotypic file for offspring:",filename,"and will store  it in population$offspring$genotypes$read\n")
		physical <-	read.table(filename,sep=sep,row.names=1,header=FALSE)
		physical <- as.matrix(physical)
		population <- intoPopulation(population, physical, "maps$physical")
		doCleanUp.internal()
	}else{
		if(verbose)cat("No physical map file:",filename,".\n")
	}
	
	#**********FINALIZING FUNCTION*************
	e <- proc.time()
	if(verbose) cat("readFiles done in",(e-s)[3],"seconds.\n")
	population$parameters$readFiles <- list(offspring, founders,sep=sep,verbose,debugMode)
	names(population$parameters$readFiles) <- c("offspring", "founders", "sep", "verbose", "debugMode")
	class(population) <- "population"
	doCleanUp.internal()
	invisible(population)
}

############################################################################################################
#									*** mapMarkers.internal ***
#
# DESCRIPTION:
#	removes from matrix1 cols or rows, which are not present in second (coparing using col/rownames)
# 
# PARAMETERS:
# 	expressionMatrix1, expressionMatrix2 - matrices with data of any type
# 	mapMode - 1 - map rows, 2 - map cols
# 	verbose - Be verbose
# 	debugMode - 1: Print our checks, 2: print additional time information
#
# OUTPUT:
#	object of class population 
#
############################################################################################################
mapMarkers.internal <- function(expressionMatrix1, expressionMatrix2, mapMode=2, verbose=FALSE, debugMode=0){
	if(mapMode==1) {
		nrRows <- nrow(expressionMatrix1)
		### warnings when names are mismatching
		if(verbose && debugMode==2)if(nrRows!=nrow(expressionMatrix2)){
			cat("Following markers will be removed:\n")
			cat(paste(rownames(expressionMatrix1)[which(!(rownames(expressionMatrix1) %in% rownames(expressionMatrix2)))],"\n"))
		}
		### mapping itself
		expressionMatrix1 <- expressionMatrix1[which(rownames(expressionMatrix1) %in% rownames(expressionMatrix2)),]
		if(verbose) cat("Because of names mismatch,",nrRows-nrow(expressionMatrix1),"markers were removed, run function with verbose=T debugMode=2 to print their names out.\n")
	}
	else if(mapMode==2){
		nrCols <- ncol(expressionMatrix1)
		### warnings when names are mismatching
		if(verbose && debugMode==2)if(nrCols!=ncol(expressionMatrix2)){
			cat("Following individuals will be removed:\n")
			paste(colnames(expressionMatrix1)[which(!(colnames(expressionMatrix1) %in% colnames(expressionMatrix2)))],"\n")
		}
		### mapping itself
		expressionMatrix1 <- expressionMatrix1[,which(colnames(expressionMatrix1) %in% colnames(expressionMatrix2))]
		if(verbose) cat("Because of names mismatch,",nrCols-ncol(expressionMatrix1),"individuals were removed, run function with verbose=T debugMode=2 to print their names out.\n")
	}
	invisible(expressionMatrix1)
}
