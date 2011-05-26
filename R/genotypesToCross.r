#####################################################################
#
# genotypesToCross.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified April 2011
# last modified in version: 0.4.3
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
#				writePhenotypes.internal, writeGenotypes.internal
#
#####################################################################



############################################################################################################
#genotypesToCross.internal- produces from genotypic matrix file containing object of type cross, reads it 
# into R a returns
# 
# ril - Ril type object, must contain parental phenotypic data.
# use - save "real" gentypes, "simulated" genotypes otr simulated genotypes ordered using "map" from gff file
# outputFile - file where object of type cross is being saved
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
#### Not yet in use:
# #doClustering - this three to be used in the future
# #groups 
# #iterations 
#
#### To be removed:
# #use=map - now map is saved in every cross, so we can use it
#
#### No longer in use:
# #limit - how many phenotypes should be written
#
############################################################################################################
#genotypesToCross.internal <- function(ril, use=c("real","simulated","map"), limit=10, doClustering=FALSE, groups=10, iterations = 100, outputFile="mycross.csv", verbose=FALSE, debugMode=0){
genotypesToCross.internal <- function(ril, use=c("real","simulated","map"), outputFile="mycross.csv", verbose=FALSE, debugMode=0){
	###CHECKS
	if(verbose && debugMode==1) cat("genotypesToCross starting.\n")
	s <- proc.time()
	
	
#**********WRITING PHENOTYPIC DATA TO FILE*************
	if(!is.null(ril$rils$phenotypes)){
		#there is phenotypic matrix
		cat("Writing phenotypic data to cross file\n")
		writePhenotypes.internal(ril, use, outputFile, verbose, debugMode)
	}else{
		#there is no phenotypic matrix
		stop("genotypesToCross not provided with phenotypic matrix, stopping\n")
	}
	
#**********WRITING GENOTYPIC DATA TO FILE*************
	if(use=="real"){
		if(is.null(ril$rils$genotypes$read)){
			stop("Use = real chosen, but there is no real genotypic data in ril$rils$genotypes$read\n")
		}else{
			cat("Cross object will be written using real genotypic data\n")
			writeGenotypes.internal(ril$rils$genotypes$read, chr=1,outputFile=outputFile, verbose=verbose, debugMode=debugMode)
		}
	}
	else if(use=="simulated" || use=="map"){
		if(is.null(ril$rils$genotypes$simulated)){
			stop("Use = simulated chosen, but there is no simulated genotypic data in ril$rils$genotypes$simulated\n")
	}else if(use=="map"){
		if(is.null(ril$rils$map)){
			stop("Use = map chosen, but there is no map data in ril$rils$map\n")
		}else{
			cat("Cross object will be written using simulated genotypic data orderd by gff map\n")
			writeGenotypesMap.internal(ril, outputFile, verbose, debugMode)
		}
	}else{
			cat("Cross object will be written using simulated genotypic data\n")
			writeGenotypes.internal(ril$rils$genotypes$simulated, chr=1, positions=ril2$location[which(ril2$location[,1] %in% ril$rils$phenotypes[(rownames(ril$rils$genotypes$simulated)),1]),2], outputFile=outputFile, verbose=verbose, debugMode=debugMode)
		}
	}	

#**********READING CROSS FILE TO R*************
	cross <- invisible(read.cross("csvr",file=outputFile, genotypes=c(0,1)))
	class(cross)[1] <- "riself"
	e <- proc.time()
	if(verbose) cat("genotypesToCross done in",(e-s)[3],"seconds.\n")
	invisible(cross)
}

############################################################################################################
#writePhenotypes.internal - sub function of genotypesToCross - writes phenotypes to file
# 
# ril - Ril type object, must contain parental phenotypic data.
# use - save "real" gentypes, "simulated" genotypes otr simulated genotypes ordered using "map" from gff file
# outputFile - file where object of type cross is being saved
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
# REMOVED: expressionMatrix - columns: individuals, rows: markers;
# 		   limit - how many phenotypes should be written
#
############################################################################################################
writePhenotypes.internal <- function(ril, use, outputFile, verbose=FALSE, debugMode=0){
	sl <- proc.time()
	
	if(use=="real"){
		if(is.null(ril$rils$genotypes$read)){
			stop("Use = real chosen, but there is no real genotypic data in ril$rils$genotypes$read\n")
		}else{
			ril$rils$phenotypes <- mapMarkers.internal(ril$rils$phenotypes,ril$rils$genotypes$read,mapMode=2)
			ril$rils$genotypes$read <- mapMarkers.internal(ril$rils$genotypes$read,ril$rils$phenotypes,mapMode=2)
			chosenLabels <- rownames(ril$rils$genotypes$read)
		}
	}else if(use=="simulated" || use=="map"){
		if(is.null(ril$rils$genotypes$simulated)){
			stop("Use = simulated or map chosen, but there is no simulated genotypic data in ril$rils$genotypes$simulated\n")
		}else{
			ril$rils$phenotypes <- mapMarkers.internal(ril$rils$phenotypes,ril$rils$genotypes$simulated,mapMode=2)
			ril$rils$genotypes$simulated <- mapMarkers.internal(ril$rils$genotypes$simulated,ril$rils$phenotypes,mapMode=2)
		}
		
		if(use=="simulated"){
				chosenLabels <- rownames(ril$rils$genotypes$simulated)
		}else if(use=="map"){
			if(is.null(ril$rils$genotypes$simulated)||is.null(ril$rils$map)){
				stop("Use =  map chosen, but there is no simulated genotypic data or map data\n")
			}else{
				chosenLabels <- rownames(ril$rils$genotypes$simulated[which(rownames(ril$rils$genotypes$simulated[,]) %in% rownames(ril$rils$map)),])
			}
		}
		
	}
	
	toWrite <- ril$rils$phenotypes[chosenLabels,]
	toWrite <- cleanNames.internal(toWrite)
	if(verbose && debugMode==1) cat("writePhenotypes starting.\n")
	write.table(cbind("","",toWrite),file=outputFile,sep=",",quote=FALSE,col.names=FALSE)
	el <- proc.time()
	if(verbose && debugMode==2)cat("Writing phenotypes done in:",(el-sl)[3],"seconds.\n")
}

############################################################################################################
#writeGenotypesMap.internal - sub function of genotypesToCross - orders markers using map and writes them to
# file using writeGenotypes.internal 
# 
# ril - Ril type object, must contain parental phenotypic data.
# outputFile - file where object of type cross is being saved
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
writeGenotypesMap.internal <- function(ril, outputFile, verbose=FALSE, debugMode=0){
	#TODO: use map to sort out the genotype, then also comparing reco map with physical
	sl <- proc.time()
	if(verbose && debugMode==1) cat("writeGenotypesMap.internal starting.\n")
	for(i in 1:length(table(ril$rils$map[,1]))){
		selectedMap <- ril$rils$map[which(ril$rils$map[,1]==i),]
		selectedGenes <- ril$rils$genotypes$simulated[which(rownames(ril$rils$genotypes$simulated) %in% rownames(selectedMap)),]
		selectedMap <- selectedMap[which(rownames(ril$rils$genotypes$simulated) %in% rownames(selectedMap)),]
		writeGenotypes.internal (selectedGenes,i,selectedMap[,2],outputFile,verbose,debugMode)
	}
	el <- proc.time()
	if(verbose && debugMode==2) cat("Writing genotypes done in:",(el-sl)[3],"seconds.\n")
}

############################################################################################################
#writeGenotypes.internal - sub function of genotypesToCross and writeGenotypesMap - writes genotypes 
# (one chromosome at the time) to file
# 
# genotypeMatrix - matrix of genotypic data, rows - markers, cols - individuals
# chr - chromosome currently being written
# outputFile - file where object of type cross is being saved
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
writeGenotypes.internal <- function(genotypeMatrix,chr=1,positions=NULL,outputFile,verbose=FALSE,debugMode=0){
	sl <- proc.time()
	cat(positions,"\n")
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
#cleanNames.internal - changing names that will crush read.cross
# 
# matrixToBeCleaned - matrix of any data type
#
############################################################################################################
cleanNames.internal <-function(matrixToBeCleaned){
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
