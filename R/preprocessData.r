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
# last modified in version: 0.6.1
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
# ril - Ril type object, must contain parental phenotypic data.
# groupLabels - Specify which column of parental data belongs to group 0 and which to group 1.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
preprocessData <- function(ril,groupLabels=c(0,0,1,1),verbose=FALSE,debugMode=0,...){
	s2<-proc.time()
	result <- NULL
	filename <- "rilrankProdRes"
	nr <- 1
	loaded <- 0
	while(file.exists(paste(filename,".RData",sep=""))&&loaded==0){
		print(paste(filename,".RData",sep=""))
		if(verbose) cat("File rilrankProdRes.Rdata already exists, reading it.\n")
		load(paste(filename,".RData",sep=""))
		if(length(result)==3){
			if(checkDataStamp.internal(ril,result[[3]])){
				ril$parental$RP <- result[[1]]
				if(verbose) cat("File rilrankProdRes.Rdata contains groupLabels used during processing. They are following:",result[[2]]," and will be used instead of ones supplied by user:",groupLabels,"\n")
				ril$parental$groups <- result[[2]]
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
		ril <- conductRP.internal(ril, groupLabels, filename, ...)
	}	
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Data preprocessing done in:",(e2-s2)[3],"seconds.\n")
	ril$parameters$preprocessData <- list("object of ril class", groupLabels,verbose,debugMode)
	names(ril$parameters$preprocessData) <- c("ril", "groupLabels", "verbose", "debugMode")
	class(ril) <- "ril"
	invisible(ril)
}

############################################################################################################
#preprocessData: Using Rank Product analysis to select differentially expressed genes.
# 
# ril - Ril type object, must contain parental phenotypic data.
# groupLabels - Specify which column of parental data belongs to group 0 and which to group 1.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
conductRP.internal <- function(ril, groupLabels, filename, ...){
	rankProdRes <- invisible(RP(ril$parental$phenotypes,groupLabels,...))
	result <- list(rankProdRes,groupLabels,makeDataStamp.internal(ril))
	save(file=paste(filename,".RData",sep=""),result)
	ril$parental$RP <- rankProdRes
	ril$parental$groups <- groupLabels
	invisible(ril)
}

############################################################################################################
#preprocessData: Using Rank Product analysis to select differentially expressed genes.
# 
# ril - Ril type object, must contain parental phenotypic data.
# groupLabels - Specify which column of parental data belongs to group 0 and which to group 1.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
makeDataStamp.internal <- function(ril){
	result <- NULL
	result$wd <- getwd()
	result$parental <- ril$parameters$readFiles$parental
	result$children <- ril$parameters$readFiles$rils
	invisible(result)
}
