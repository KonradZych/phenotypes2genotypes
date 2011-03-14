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

genotypesToCross <- function(genotypeMatrix, expressionMatrix, doClustering=FALSE, groups=10, iterations = 100, outputFile="mycross.csv", verbose=FALSE, debugMode=0){
	###CHECKS
	
	if(verbose && debugMode==1) cat("genotypesToCross starting.\n")
	s <- proc.time()
	cat("",file=outputFile)
	#saving phenotypic data
	writePhenotypes(expressionMatrix, outputFile, verbose, debugMode)
	if(!doClustering){
		groups <- 1
	}
	if(groups>1){
		r <- bestClustering(genotypeMatrix,groups,iterations,verbose,debugMode)
		r <- kmeans(r,groups)
	}else{
		r <- kmeans(genotypeMatrix,groups)
	}
	sorted <- sort(r[[1]],index.return=TRUE)
	for(i in 1:groups){
		sl <- proc.time()
		if(verbose && debugMode==1) cat("writeGenotypes starting  for chromome",i,"out of",groups,".\n")
		writeGenotypes(genotypeMatrix[sorted[[2]][which(sorted[[1]]==i)],], i,outputFile, verbose, debugMode)
		el <- proc.time()
		if(verbose && debugMode==2)cat("writeGenotypes for chromome",i," done in:",(el-sl)[3],"seconds.\n")
	}  
	#reading freshly made file to R
	cross <- invisible(read.cross("csvr",file=outputFile, genotypes=c(0,1)))
	#forcing cross time to RIL
	class(cross)[1] <- "riself"
	e <- proc.time()
	#returning cross
	if(verbose) cat("genotypesToCross done in",(e-s)[3],"seconds.\n")
	invisible(cross)
}

#writePhenotypes - writes to file phenotypic data (cross object format)
writePhenotypes <- function(expressionMatrix, outputFile, verbose=FALSE, debugMode=TRUE){
	sl <- proc.time()
	if(verbose && debugMode==1) cat("writePhenotypes starting.\n")
	write.table(cbind("","",expressionMatrix),file=outputFile,sep=",",quote=FALSE,col.names=FALSE)
	el <- proc.time()
	if(verbose && debugMode==2)cat("Writing genotypes done in:",(el-sl)[3],"seconds.\n")
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
