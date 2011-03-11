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
# usage cross <- orderedCross(genotypicMatrix)

genotypesToCross <- function(genotypeMatrix, expressionMatrix, doClustering=FALSE, groups=10, outputFile="mycross.csv", verbose=FALSE, debugMode=0){
	if(verbose) cat("genotypesToCross starting.\n\n")
	s <- proc.time()
	
	#printing faked phenotype
	cat("",file=outputFile)
	cat("phenotype",sep="",file=outputFile,append=TRUE)
	cat(",,,",sep="",file=outputFile,append=TRUE)
	cat(paste(runif(48,0,100),collapse=","),sep="",file=outputFile,append=TRUE)
	cat("\n",sep="",file=outputFile,append=TRUE)
	
	#Dividing data in ten groups based on corelation using kmeans
	cor_matrix <- cor(t(genotypeMatrix), use="pairwise.complete.obs")
	r <- kmeans(genotypeMatrix,groups)
	sorted <- sort(r[[1]],index.return=TRUE)

	#printing genotypic data to file in specified format
	for(i in 1:groups){
		sl <- proc.time()
		if(verbose){
			cat("   ### Writing chromosome: ",i,"nr of markers:",length(which(r[[1]]==i)),"###\n")
			cat(paste("      -> Marker ",colnames(sorted[[1]][which(sorted[[1]]==i)]),"\n",sep=""))
		}
		write.table(cbind(sorted[[1]][which(sorted[[1]]==i)], 1:table(sort(r[[1]]))[i], genotypeMatrix[sorted[[2]][which(sorted[[1]]==i)],]) , file=outputFile , sep="," , quote=FALSE , col.names=FALSE , append=TRUE)
		el <- proc.time()
		if(verbose){cat("   ### Writing chromosome:",i,"taken:",(el-sl)[3],"seconds. ###\n\n\n")}
	}
  
	#reading freshly made file to R
	cross <- read.cross("csvr",file=outputFile,genotypes=c(0,1))
	#forcing cross time to RIL
	class(cross)[1] <- "riself"
	e <- proc.time()
	if(verbose){cat("Done without errors in:",(e-s)[3],"seconds.\n")}
	#returning cross
	invisible(cross)
}



test.genotypesToCross <- function(){

}
