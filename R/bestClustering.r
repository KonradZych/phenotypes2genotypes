#####################################################################
#
# bestClustering.R
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

#bestClustering
bestClustering <- function(genotypeMatrix, groups, iterations, verbose=FALSE, debugMode=0){
	#CHECKS
	s <- proc.time()
	if(verbose && debugMode==1) cat("genotypesToCross starting.\n")
	
	#clustering
	recoMatrix <- recombinationCount(genotypeMatrix,flip=1)
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
	pointsMatrix <- apply(groupsMatrix,2,bestClusteringSub,groupsMatrix)
	ep <- proc.time()
	if(verbose && debugMode==2)cat("Pointing done in:",(ep-sp)[3],"seconds.\n")
	
	e <- proc.time()
	if(verbose) cat("bestClustering done in",(e-s)[3],"seconds.\n")
	invisible(pointsMatrix)
}


bestClusteringSub <- function(groupsMatrixCol, groupsMatrix, verbose=FALSE, debugMode=0){
	#CHECKS
	s <- proc.time()
	if(verbose && debugMode==1) cat("genotypesToCross starting.\n")
	pointsMatrixCol <- apply(groupsMatrix,2,bestClusteringSubSub,groupsMatrixCol)
	e <- proc.time()
	if(verbose) cat("bestClusteringSub done in",(e-s)[3],"seconds.\n")
	invisible(pointsMatrixCol)
}

bestClusteringSubSub <- function(groupsMatrixCol1,groupsMatrixCol2){
	point <- sum(groupsMatrixCol1==groupsMatrixCol2)
	invisible(point)
}