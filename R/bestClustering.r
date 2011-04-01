#####################################################################
#
# bestClustering.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified April 2011
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
# Contains: bestClustering
#				bestClusteringSub, bestClusteringSubSub
#
#####################################################################

############################################################################################################
#bestClustering- performs clustering using kmeans in specified iteriations number.
# 
# genotypeMatrix - Matrix with genotype values with: columns, individuals and on the rows, markers
# groups - number of groups data should be splitted to
# iterations - ow many times clustering should be performed
# use - r -> use recombination matrix for clustering, c -> corelation matrix
# flip - specifies whether one of the rows that are being compared should be flipped(1) or not(0)
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
bestClustering.internal <- function(genotypeMatrix, groups, iterations, use="r", flip=0, verbose=FALSE, debugMode=0){
	#CHECKS
	s <- proc.time()
	if(verbose && debugMode==1) cat("bestClustering starting.\n")
	
	#clustering
	if(use=="r"){
		sr <- proc.time()
		recoMatrix <- recombinationCount(genotypeMatrix,flip=flip)
		er <- proc.time()
		if(verbose && debugMode==2)cat("Recombination matrix created, done in:",(er-sr)[3],"seconds.\n")
	}else if(use=="c"){
		sr <- proc.time()
		recoMatrix <- cor(genotypeMatrix,use="pairwise.complete.obs")
		er <- proc.time()
		if(verbose && debugMode==2)cat("Corelation matrix created, done in:",(er-sr)[3],"seconds.\n")
	}
	groupsMatrix <- NULL
	si <- proc.time()
	for(i in 1:iterations){
		#ngroups <- groups + sample(3,1)-2
		currentGroups <- kmeans(recoMatrix,groups)
		groupsMatrix <- rbind(groupsMatrix,(as.numeric(currentGroups[[1]])))
	}
	ei <- proc.time()
	if(verbose && debugMode==2)cat("Clustering using kmeans done",iterations," times, done in:",(ei-si)[3],"seconds.\n")
	
	#marking
	sp <- proc.time()
	pointsMatrix <- apply(groupsMatrix,2,bestClusteringSub.internal,groupsMatrix)
	ep <- proc.time()
	if(verbose && debugMode==2)cat("Pointing done in:",(ep-sp)[3],"seconds.\n")
	
	e <- proc.time()
	if(verbose) cat("bestClustering done in",(e-s)[3],"seconds.\n")
	invisible(pointsMatrix)
}

############################################################################################################
#bestClusteringSub - sub function of bestClustering - comparing single column vs all the others using 
# bestClusteringSubSub
# 
# genotypeMatrixCol - Vector with genotype values for single individual
# groupsMatrix - Matrix cols: individuals, rows: clustering iterations 
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
bestClusteringSub.internal <- function(groupsMatrixCol, groupsMatrix, verbose=FALSE, debugMode=0){
	#CHECKS
	s <- proc.time()
	if(verbose && debugMode==1) cat("genotypesToCross starting.\n")
	pointsMatrixCol <- apply(groupsMatrix,2,bestClusteringSubSub.internal,groupsMatrixCol)
	e <- proc.time()
	if(verbose) cat("bestClusteringSub done in",(e-s)[3],"seconds.\n")
	invisible(pointsMatrixCol)
}

############################################################################################################
#bestClusteringSubSub - comparing two vectors
# 
# genotypeMatrixCol1 - Vector with genotype values for single individual
# genotypeMatrixCol2 - Vector with genotype values for single individual
#
############################################################################################################
bestClusteringSubSub.internal <- function(groupsMatrixCol1,groupsMatrixCol2){
	point <- sum(groupsMatrixCol1==groupsMatrixCol2)
	invisible(point)
}