#################################################################################
#
# preprocessData.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified June 2011
# last modified in version: 0.7.1
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
# Contains: preprocessData 
#
#################################################################################

############################################################################################################
#preprocessData: Using Rank Product analysis to select differentially expressed genes.
# 
# population - Ril type object, must contain founders phenotypic data.
# groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
preprocessData <- function(population,groupLabels=c(0,0,1,1),verbose=FALSE,debugMode=0,...){
	s2<-proc.time()
	result <- NULL
	filename <- "populationRankProdRes"
	nr <- 1
	loaded <- 0
	while(file.exists(paste(filename,".RData",sep=""))&&loaded==0){
		print(paste(filename,".RData",sep=""))
		if(verbose) cat("File populationRankProdRes.Rdata already exists, reading it.\n")
		load(paste(filename,".RData",sep=""))
		if(length(result)==3){
			if(checkDataStamp.internal(population,result[[3]])){
				population$founders$RP <- result[[1]]
				if(verbose) cat("File populationRankProdRes.Rdata contains groupLabels used during processing. They are following:",result[[2]]," and will be used instead of ones supplied by user:",groupLabels,"\n")
				population$founders$groups <- result[[2]]
				loaded <- 1
			}else{
			filename <- paste(filename,nr,sep="_")
			nr <- nr+1
			}
		}else{
			filename <- paste(filename,nr,sep="_")
			nr <- nr+1
		}
	}
	if(loaded==0){
		population <- conductRP.internal(population, groupLabels, filename, ...)
	}	
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Data preprocessing done in:",(e2-s2)[3],"seconds.\n")
	population$parameters$preprocessData <- list("object of population class", groupLabels,verbose,debugMode)
	names(population$parameters$preprocessData) <- c("population", "groupLabels", "verbose", "debugMode")
	class(population) <- "population"
	invisible(population)
}

############################################################################################################
#preprocessData: Using Rank Product analysis to select differentially expressed genes.
# 
# population - Ril type object, must contain founders phenotypic data.
# groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
conductRP.internal <- function(population, groupLabels, filename, ...){
	rankProdRes <- invisible(RP(population$founders$phenotypes,groupLabels,...))
	result <- list(rankProdRes,groupLabels,makeDataStamp.internal(population))
	save(file=paste(filename,".RData",sep=""),result)
	population$founders$RP <- rankProdRes
	population$founders$groups <- groupLabels
	invisible(population)
}

############################################################################################################
#preprocessData: Using Rank Product analysis to select differentially expressed genes.
# 
# population - Ril type object, must contain founders phenotypic data.
# groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
makeDataStamp.internal <- function(population){
	result <- NULL
	result$wd <- getwd()
	if(nrow(population$founders$phenotypes)<10){
		rowS <- nrow(population$founders$phenotypes)
	}else{
		rowS <- 10
	}
	result$rowS <- rowS
	result$founders <- population$founders$phenotypes[1:rowS,]
	invisible(result)
}
