#############################################################################################
#
# genotypesToCross.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified May 2011
# last modified in version: 0.7.1
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
#genotypesToCross.internal- produces from genotypic matrix file containing object of type cross, reads it 
# into R a returns
# 
# population - population type object, must contain founders phenotypic data.
# use - save "real" gentypes, "simulated" genotypes otr simulated genotypes ordered using "map" from gff file
# outputFile - file where object of type cross is being saved
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
#
############################################################################################################
genotypesToCross.internal <- function(population, use=c("simulated","map_genetic","map_physical","real"), outputFile="mycross.csv", verbose=FALSE, debugMode=0){
	###CHECKS
	if(verbose && debugMode==1) cat("genotypesToCross starting.\n")
	s <- proc.time()
	
	
#**********WRITING PHENOTYPIC DATA TO FILE*************
	if(!is.null(population$offspring$phenotypes)){
		#there is phenotypic matrix
		cat("Writing phenotypic data to cross file\n")
		writePhenotypes.internal(population, use, outputFile, verbose, debugMode)
	}else{
		#there is no phenotypic matrix
		stop("genotypesToCross not provided with phenotypic matrix, stopping\n")
	}
	
#**********WRITING GENOTYPIC DATA TO FILE*************
	if(use=="real"){
		if(is.null(population$offspring$genotypes$real)){
			stop("Use = real chosen, but there is no real genotypic data in population$offspring$genotypes$read\n")
		}else{
			cat("Cross object will be written using real genotypic data\n")
			writeGenotypes.internal(population$offspring$genotypes$read, chr=1,outputFile=outputFile, verbose=verbose, debugMode=debugMode)
		}
	}
	else if(use=="simulated" || use=="map_genetic" || use=="map_physical"){
		if(is.null(population$offspring$genotypes$simulated)){
				stop("Use = simulated chosen, but there is no simulated genotypic data in population$offspring$genotypes$simulated\n")
		}else if(use=="map_genetic"){
			### check this, doesn't look nice!
			if(is.null(population$maps$genetic)){
				stop("Use = map chosen, but there is no map data in population$offspring$map\n")
			}else{
				cat("Cross object will be written using simulated genotypic data orderd by gff map\n")
				writeGenotypes.internal(population$offspring$genotypes$simulated, chr=(population$maps$genetic[rownames(population$offspring$genotypes$simulated),1], positions=(population$maps$genetic[rownames(population$offspring$genotypes$simulated),2], outputFile=outputFile, verbose=verbose, debugMode=debugMode)
			}
		}else if(use=="map_physical"){
			### check this, doesn't look nice!
			if(is.null(population$maps$physical)){
				stop("Use = map chosen, but there is no map data in population$offspring$map\n")
			}else{
				cat("Cross object will be written using simulated genotypic data orderd by gff map\n")
				writeGenotypes.internal(population$offspring$genotypes$simulated, chr=(population$maps$physical[rownames(population$offspring$genotypes$simulated),1], positions=(population$maps$physical[rownames(population$offspring$genotypes$simulated),2], outputFile=outputFile, verbose=verbose, debugMode=debugMode)
			}
		}else{
				cat("Cross object will be written using simulated genotypic data\n")
				writeGenotypes.internal(population$offspring$genotypes$simulated, chr=1, outputFile=outputFile, verbose=verbose, debugMode=debugMode)
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
# population - Ril type object, must contain founders phenotypic data.
# use - save "real" gentypes, "simulated" genotypes otr simulated genotypes ordered using "map" from gff file
# outputFile - file where object of type cross is being saved
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
# REMOVED: expressionMatrix - columns: individuals, rows: markers;
# 		   limit - how many phenotypes should be written
#
############################################################################################################
writePhenotypes.internal <- function(population, use, outputFile, verbose=FALSE, debugMode=0){
	sl <- proc.time()
	
	if(use=="real"){
		if(is.null(population$offspring$genotypes$read)){
			stop("Use = real chosen, but there is no real genotypic data in population$offspring$genotypes$read\n")
		}else{
			population$offspring$phenotypes <- mapMarkers.internal(population$offspring$phenotypes,population$offspring$genotypes$read,mapMode=2)
			population$offspring$genotypes$read <- mapMarkers.internal(population$offspring$genotypes$read,population$offspring$phenotypes,mapMode=2)
			chosenLabels <- rownames(population$offspring$genotypes$read)
		}
	}else if(use=="simulated" || use=="map"){
		if(is.null(population$offspring$genotypes$simulated)){
			stop("Use = simulated or map chosen, but there is no simulated genotypic data in population$offspring$genotypes$simulated\n")
		}else{
			population$offspring$phenotypes <- mapMarkers.internal(population$offspring$phenotypes,population$offspring$genotypes$simulated,mapMode=2)
			population$offspring$genotypes$simulated <- mapMarkers.internal(population$offspring$genotypes$simulated,population$offspring$phenotypes,mapMode=2)
		}
		
		if(use=="simulated"){
				chosenLabels <- rownames(population$offspring$genotypes$simulated)
		}else if(use=="map"){
			if(is.null(population$offspring$genotypes$simulated)||is.null(population$offspring$map)){
				stop("Use =  map chosen, but there is no simulated genotypic data or map data\n")
			}else{
				chosenLabels <- rownames(population$offspring$genotypes$simulated[which(rownames(population$offspring$genotypes$simulated[,]) %in% rownames(population$offspring$map)),])
			}
		}
		
	}
	
	toWrite <- population$offspring$phenotypes[chosenLabels,]
	toWrite <- cleanNames.internal(toWrite)
	if(verbose && debugMode==1) cat("writePhenotypes starting.\n")
	write.table(cbind("","",toWrite),file=outputFile,sep=",",quote=FALSE,col.names=FALSE)
	el <- proc.time()
	if(verbose && debugMode==2)cat("Writing phenotypes done in:",(el-sl)[3],"seconds.\n")
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
