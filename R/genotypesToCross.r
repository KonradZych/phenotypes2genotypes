#####################################################################
#
# genotypesToCross.R
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
# Contains: genotypesToCross
#				writePhenotypes, doClustering, writeGenotypes, fakeGenotypes, fakePhenotypes
#
#####################################################################


#genotypesToCross - produces from genotypic matrix file containing object of type cross, reads it into R a returns
# genotypeMatrix - matrix of genotypic data, rows - markers, cols - individuals
# expressionMatrix - columns: individuals, rows: markers
# nr_iterations - not yet used, maybe never will be
# groups - nr of groups we are dividing our dataset to, hopefully coresponds to nr of chromosomes
# outputFile - file where object of type cross is being saved
# verbose - standard
# debugMode - standard
# genos - argument passed to read.cross (chars describing genotypes)
# usage cross <- orderedCross(genotypicMatrix,expressionMatrix)

genotypesToCross <- function(ril, doClustering=FALSE, groups=10, iterations = 100, outputFile="mycross.csv", verbose=FALSE, debugMode=0){
	###CHECKS
	if(verbose && debugMode==1) cat("genotypesToCross starting.\n")
	s <- proc.time()
	
	#**********CHECKING GENOTYPIC MATRICES*************
	cat("",file=outputFile)
	if(!is.null(ril$rils$genotypes$read)){
		if(!is.null(ril$rils$genotypes$simulated)){
			#two matrices present
			crossMode <- 0
		}
		#only read from file matrix present
		crossMode <- 1
	}else if(!is.null(ril$rils$genotypes$simulated)){
		#only simulated matrix presen
		crossMode <- 2
	}else{
		cat("genotypesToCross not provided with genotypic matrix, will fake one!\n")
		fakeGenotypes(outputFile, verbose, debugMode)
		crossMode <- 3
	}
	
	#**********WRITING PHENOTYPIC DATA TO FILE*************
	if(!is.null(ril$rils$phenotypes)){
		#there is phenotypic matrix
		writePhenotypes(expressionMatrix, outputFile, verbose, debugMode)
	}else if(crossMode!=3){
		#there is no phenotypic matrix, but there is genotypic matrix
		fakePhenotypes(outputFile, verbose, debugMode)
	}else{
		#there is no phenotypic matrix, neither the genotypic matrix
		stop("genotypesToCross not provided with genotypic neither phenotypic matrix, stopping\n")
	}

	#**********WRITING GENOTYPIC DATA TO FILE*************
	if(doClustering){
		if(crossMode==0){
			doClustering(ril$rils$genotypes$read,groups,outputFile, verbose, debugMode)
			doClustering(ril$rils$genotypes$simulated,groups,outputFile, verbose, debugMode)
		}else if(crossMode==1){
			doClustering(ril$rils$genotypes$read,groups,outputFile, verbose, debugMode)
		}else if(crossMode==2){
			doClustering(ril$rils$genotypes$simulated,groups,outputFile, verbose, debugMode)
		}
	}else{
		if(crossMode==0){
			writeGenotypes(ril$rils$genotypes$read, 1,outputFile, verbose, debugMode)
			writeGenotypes(ril$rils$genotypes$simulated, 1,outputFile, verbose, debugMode)
		}else if(crossMode==1){
			writeGenotypes(ril$rils$genotypes$read, 1,outputFile, verbose, debugMode)
		}else if(crossMode==2){
			writeGenotypes(ril$rils$genotypes$simulated, 1,outputFile, verbose, debugMode)
		}
	}

	#**********READING CROSS FILE TO R*************
	cross <- invisible(read.cross("csvr",file=outputFile, genotypes=c(0,1)))
	class(cross)[1] <- "riself"
	e <- proc.time()
	if(verbose) cat("genotypesToCross done in",(e-s)[3],"seconds.\n")
	invisible(cross)
}

#writePhenotypes - writes to file phenotypic data (cross object format)
writePhenotypes <- function(expressionMatrix, outputFile, verbose=FALSE, debugMode=TRUE){
	sl <- proc.time()
	if(verbose && debugMode==1) cat("writePhenotypes starting.\n")
	write.table(cbind("","",expressionMatrix),file=outputFile,sep=",",quote=FALSE,col.names=FALSE)
	el <- proc.time()
	if(verbose && debugMode==2)cat("Writing phenotypes done in:",(el-sl)[3],"seconds.\n")
}

doClustering <- function(genotypeMatrix, groups, iterations, outputFile, verbose=FALSE, debugMode=0){
	r <- bestClustering(genotypeMatrix,groups,iterations,verbose=verbose,debugMode=debugMode)
	r <- kmeans(r,groups)
	sorted <- sort(r[[1]],index.return=TRUE)
	for(i in 1:groups){
		sl <- proc.time()
		if(verbose && debugMode==1) cat("writeGenotypes starting  for chromome",i,"out of",groups,".\n")
		writeGenotypes(genotypeMatrix[sorted[[2]][which(sorted[[1]]==i)],], i,outputFile, verbose, debugMode)
		el <- proc.time()
	}
	if(verbose && debugMode==2)cat("doClustering done in:",(el-sl)[3],"seconds.\n")
}


#writeGenotypes - writes to file genotypic data (cross object format)
#genotypeMatrix - matrix of genotypic data, rows - markers, cols - individuals
#verbose - standard
#debugMode - standard 1 -> gives info, that function is starting  2 -> gives additional time information
writeGenotypes <- function(genotypeMatrix,chr=1,outputFile,verbose=FALSE,debugMode=0){
	sl <- proc.time()
	if(verbose && debugMode==1) cat("writeGenotypes starting.\n")
	write.table(cbind(chr,1:nrow(genotypeMatrix),genotypeMatrix),file=outputFile,sep=",",quote=FALSE,col.names=FALSE,append=TRUE)
	el <- proc.time()
	if(verbose && debugMode==2) cat("Writing genotypes done in:",(el-sl)[3],"seconds.\n")
}

fakeGenotypes <- function(len, outputFile, verbose=FALSE, debugMode=0){
	cat("faked,1,1,",file=outputFile,append=TRUE)
	cat(paste(sample(c(0,1),100,T),",",sep="",collapse=""),file=outputFile,append=TRUE)
	cat("\n",file=outputFile,append=TRUE)
}

fakePhenotypes <- function(len, outputFile, verbose=FALSE, debugMode=0){
	cat("faked,,,",file=outputFile,append=TRUE)
	cat(paste(sample(c(0.03,0.04,0.05),100,T),",",sep="",collapse=""),file=outputFile,append=TRUE)
	cat("\n",file=outputFile,append=TRUE)
}
