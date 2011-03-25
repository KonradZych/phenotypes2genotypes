#################################################################################
#
# preprocessData.R
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
# Contains: preprocessData 
# 				mapMarkers
#
#################################################################################

#preprocessData
preprocessData <- function(ril,groupLabels=c(0,0,1,1),verbose=FALSE,debugMode=0,...){
	s2<-proc.time()
	require(RankProd)
	if(file.exists("rilRP.Rdata")){
		if(verbose) cat("File rilRP.Rdata already exists, reading it.\n")
		load("rilRP.Rdata")
		ril$parental$RP <- res
	}else{
		#wasting memory here because of Rbug
		res <- invisible(RP(ril$parental$phenotypes,groupLabels,...))
		save(file="rilRP.Rdata",res)
		ril$parental$RP <- res
	}
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Data preprocessing done in:",(e2-s2)[3],"seconds.\n")
	invisible(ril)
}

mapMarkers.internal <- function(expressionMatrix1, expressionMatrix2, mapMode=2,verbose=FALSE,debugMode=0){
	if(mapMode==1) {
		nrRows <- nrow(expressionMatrix1)
		#warnings when names are mismatching
		if(verbose && debugMode==2)if(nrRows!=nrow(expressionMatrix2)){
			cat("Following markers will be removed:\n")
			cat(paste(rownames(expressionMatrix1)[which(!(rownames(expressionMatrix1) %in% rownames(expressionMatrix2)))],"\n"))
		}
		#mapping itself
		expressionMatrix1 <- expressionMatrix1[which(rownames(expressionMatrix1) %in% rownames(expressionMatrix2)),]
		if(verbose) cat("Because of names mismatch,",nrRows-nrow(expressionMatrix1),"markers were removed, run function with verbose=T debugMode=2 to print their names out.\n")
	}
	else if(mapMode==2){
		nrCols <- ncol(expressionMatrix1)
		#warnings when names are mismatching
		if(verbose && debugMode==2)if(nrCols!=ncol(expressionMatrix2)){
			cat("Following individuals will be removed:\n")
			paste(colnames(expressionMatrix1)[which(!(colnames(expressionMatrix1) %in% colnames(expressionMatrix2)))],"\n")
		}
		#mapping itself
		expressionMatrix1 <- expressionMatrix1[,which(colnames(expressionMatrix1) %in% colnames(expressionMatrix2))]
		if(verbose) cat("Because of names mismatch,",nrCols-ncol(expressionMatrix1),"individuals were removed, run function with verbose=T debugMode=2 to print their names out.\n")
	}
	invisible(expressionMatrix1)
}
