#############################################################################################
#
# genotypesToCross.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified July 2011
# last modified in version: 0.8.1
# in current version: active, internal in main workflow
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
# Contains: genotypesToCross.internal 
#				writePhenotypes.internal, writeGenotypes.internal, cleanNames.internal
#
#############################################################################################



############################################################################################################
#									*** genotypesToCross.internal ***
#
# DESCRIPTION:
#	produces from genotypic matrix file containing object of type cross, reads it into R a returns
# 
# PARAMETERS:
# 	population - population type object, must contain founders phenotypic data.
# 	use - save "real" gentypes, "simulated" genotypes otr simulated genotypes ordered using "map" from gff file
# 	outputFile - file where object of type cross is being saved
# 	verbose - Be verbose
# 	debugMode - 1: Print our checks, 2: print additional time information
#
# OUTPUT:
#	object of class cross
#
############################################################################################################
genotypesToCross.internal <- function(population, genotype=c("simulated","real"), orderUsing=c("none","map_genetic","map_physical"), outputFile="mycross.csv", verbose=FALSE, debugMode=0){
	###CHECKS
	is.population(population)
	if(verbose && debugMode==1) cat("genotypesToCross starting.\n")
	s <- proc.time()
	if(orderUsing=="map_physical"&&is.null(population$maps$physical)) stop("orderUsing=map_physical chosen, but there is no map in population$maps$physical\n")
	if(orderUsing=="map_genetic"&&is.null(population$maps$genetic)) stop("orderUsing=map_genetic chosen, but there is no map in population$maps$genetic\n")
	
#**********WRITING PHENOTYPIC DATA TO FILE*************
	if(!is.null(population$offspring$phenotypes)){
		#there is phenotypic matrix
		cat("Writing phenotypic data to cross file\n")
		population<-writePhenotypes.internal(population, genotype, outputFile, verbose, debugMode)
	}else{
		#there is no phenotypic matrix
		stop("genotypesToCross not provided with phenotypic matrix, stopping\n")
	}
	
#**********WRITING GENOTYPIC DATA TO FILE*************
	if(genotype=="real"){
		if(is.null(population$offspring$genotypes$real)){
			stop("Use = real chosen, but there is no real genotypic data in population$offspring$genotypes$real\n")
		}else{
			genoL <- length(table(population$offspring$genotypes$real))
			if(orderUsing=="none"){
				writeGenotypes.internal(population$offspring$genotypes$real, chr=1, outputFile=outputFile, verbose=verbose, debugMode=debugMode)
			}else if(orderUsing=="map_physical"){
				population$maps$physical <- mapMarkers.internal(population$maps$physical,population$offspring$genotypes$real, mapMode=1, verbose=verbose)
				if(is.null(population$maps$physical)) stop("orderUsing = map_physical chosen, but no physical map provided in population$maps$physical\n")
				writeGenotypes.internal(population$offspring$genotypes$real, chr=population$maps$physical[rownames(population$offspring$genotypes$real),1], positions=population$maps$physical[rownames(population$offspring$genotypes$real),2], outputFile=outputFile, verbose=verbose, debugMode=debugMode)
			}else if(orderUsing=="map_genetic"){
				population$maps$genetic <- mapMarkers.internal(population$maps$genetic,population$offspring$genotypes$real, mapMode=1, verbose=verbose)
				if(is.null(population$maps$genetic)) stop("orderUsing = map_physical chosen, but no genetic map provided in population$maps$genetic\n")
				writeGenotypes.internal(population$offspring$genotypes$real, chr=population$maps$genetic[rownames(population$offspring$genotypes$real),1], positions=population$maps$genetic[rownames(population$offspring$genotypes$real),2], outputFile=outputFile, verbose=verbose, debugMode=debugMode)
			}
		
		}
	}
	else if(genotype=="simulated"){
		if(is.null(population$offspring$genotypes$simulated)){
			stop("Use = simulated chosen, but there is no simulated genotypic data in population$offspring$genotypes$simulated\n")
		}else{
			genoL <- length(table(population$offspring$genotypes$simulated))
			if(orderUsing=="none"){
				writeGenotypes.internal(population$offspring$genotypes$simulated, chr=1, outputFile=outputFile, verbose=verbose, debugMode=debugMode)
			}else if(orderUsing=="map_physical"){
				population$maps$physical <- mapMarkers.internal(population$maps$physical,population$offspring$genotypes$simulated, mapMode=1, verbose=verbose)
				if(is.null(population$maps$physical)) stop("orderUsing = map_physical chosen, but no physical map provided in population$maps$physical\n")
				writeGenotypes.internal(population$offspring$genotypes$simulated, chr=population$maps$physical[rownames(population$offspring$genotypes$simulated),1], positions=population$maps$physical[rownames(population$offspring$genotypes$simulated),2], outputFile=outputFile, verbose=verbose, debugMode=debugMode)
			}else if(orderUsing=="map_genetic"){
				population$maps$genetic <- mapMarkers.internal(population$maps$genetic,population$offspring$genotypes$simulated, mapMode=1, verbose=verbose)
				if(is.null(population$maps$genetic)) stop("orderUsing = map_physical chosen, but no genetic map provided in population$maps$genetic\n")
				writeGenotypes.internal(population$offspring$genotypes$simulated, chr=population$maps$genetic[rownames(population$offspring$genotypes$simulated),1], positions=population$maps$genetic[rownames(population$offspring$genotypes$simulated),2], outputFile=outputFile, verbose=verbose, debugMode=debugMode)
			}
		}
	}	

if(genoL==2){
	genotypes <- c(0,1)
}else if(genoL==3){
	genotypes <- c(0,1,2)
}
	
#**********READING CROSS FILE TO R*************
	cross <- invisible(read.cross("csvr",file=outputFile, genotypes=genotypes))
	class(cross)[1] <- "riself"
	e <- proc.time()
	if(verbose) cat("genotypesToCross done in",(e-s)[3],"seconds.\n")
	invisible(cross)
}

############################################################################################################
#									*** writePhenotypes.internal ***
#
# DESCRIPTION:
#	sub function of genotypesToCross - writes phenotypes to file
# 
# PARAMETERS:
# 	population - Ril type object, must contain founders phenotypic data.
# 	use - save "real" gentypes, "simulated" genotypes otr simulated genotypes ordered using "map" from gff file
# 	outputFile - file where object of type cross is being saved
# 	verbose - Be verbose
# 	debugMode - 1: Print our checks, 2: print additional time information
#
# OUTPUT:
#	none
#
############################################################################################################
writePhenotypes.internal <- function(population, genotype, outputFile, verbose=FALSE, debugMode=0){
	sl <- proc.time()
	if(verbose && debugMode==1) cat("writePhenotypes starting.\n")
	if(genotype=="real"){
		if(is.null(population$offspring$genotypes$real)){
			stop("genotype = real chosen, but there is no real genotypic data in population$offspring$genotypes$read\n")
		}else{
			population$offspring$phenotypes <- mapMarkers.internal(population$offspring$phenotypes,population$offspring$genotypes$real,mapMode=2)
			population$offspring$genotypes$real <- mapMarkers.internal(population$offspring$genotypes$real,population$offspring$phenotypes,mapMode=2)
		}
	}else if(genotype=="simulated"){
		if(is.null(population$offspring$genotypes$simulated)){
			stop("genotype = simulated or map chosen, but there is no simulated genotypic data in population$offspring$genotypes$simulated\n")
		}else{
			population$offspring$phenotypes <- mapMarkers.internal(population$offspring$phenotypes,population$offspring$genotypes$simulated,mapMode=2)
			population$offspring$genotypes$simulated <- mapMarkers.internal(population$offspring$genotypes$simulated,population$offspring$phenotypes,mapMode=2)
		}
	}
	if(nrow(population$offspring$phenotypes)>1000){
		population$offspring$phenotypes <- population$offspring$phenotypes[1:1000,]
	}
	population$offspring$phenotypes<- cleanNames.internal(population$offspring$phenotypes)
	write.table(cbind("","",population$offspring$phenotypes),file=outputFile,sep=",",quote=FALSE,col.names=FALSE)
	el <- proc.time()
	if(verbose && debugMode==2)cat("Writing phenotypes done in:",(el-sl)[3],"seconds.\n")
	invisible(population)
}

############################################################################################################
#									*** writeGenotypes.internal ***
#
# DESCRIPTION:
#	sub function of genotypesToCross and writeGenotypesMap - writes genotypes (one chromosome at the time) 
#	to file
# 
# PARAMETERS:
# 	genotypeMatrix - matrix of genotypic data, rows - markers, cols - individuals
# 	chr - chromosome currently being written
# 	outputFile - file where object of type cross is being saved
# 	verbose - Be verbose
# 	debugMode - 1: Print our checks, 2: print additional time information
#
# OUTPUT:
#	none
#
############################################################################################################
writeGenotypes.internal <- function(genotypeMatrix,chr=1,positions=NULL,outputFile,verbose=FALSE,debugMode=0){
	sl <- proc.time()
	if(verbose && debugMode==1) cat("writeGenotypes starting.\n")
	if(is.null(positions)) positions <- 1:nrow(genotypeMatrix)
	else if(length(positions)!=length(1:nrow(genotypeMatrix))) stop("Posistions object is not correct, check help files.\n")
	if(verbose && debugMode==1) cat("writeGenotypes starting.\n")
	genotypeMatrix <- cleanNames.internal(genotypeMatrix)
	write.table(cbind(rownames(genotypeMatrix),chr,positions,genotypeMatrix),file=outputFile,sep=",",quote=FALSE,
		col.names=FALSE,append=TRUE,row.names=FALSE)
	el <- proc.time()
	if(verbose && debugMode==2) cat("Writing genotypes done in:",(el-sl)[3],"seconds.\n")
}

############################################################################################################
#									*** cleanNames.internal ***
#
# DESCRIPTION:
#	changing names that will crush read.cross
# 
# PARAMETERS:
# 	matrixToBeCleaned - matrix of any data type
#
# OUTPUT:
#	matrix of any data type
#
############################################################################################################
cleanNames.internal <- function(matrixToBeCleaned){
	for(i in 1:nrow(matrixToBeCleaned)){
		old <- rownames(matrixToBeCleaned)[i]
		new <- gsub(",","_",rownames(matrixToBeCleaned)[i])
		if(old != new){
			rownames(matrixToBeCleaned)[i] <- new
			cat("WARNING: marker name switched from:",old,"to",new,"because it contained ','!\n")
		}
	}
	invisible(matrixToBeCleaned)
}
