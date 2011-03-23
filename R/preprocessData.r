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
	}else{
		#wasting memory here because of Rbug
		res <- invisible(RP(ril$parental$phenotypes,groupLabels,...))
		print(ril$parental$RP$pval[[1]][2])
		save(file="rilRP.Rdata",ril$parental$RP)
		ril$parental$RP <- res
	}
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Data preprocessing done in:",(e2-s2)[3],"seconds.\n")
	invisible(ril)
}

mapMarkers <- function(expressionMatrix1, expressionMatrix2, mapMode=2){
	if(mapMode==1) expressionMatrix1 <- expressionMatrix1[which(rownames(expressionMatrix1) %in% rownames(expressionMatrix2)),]
	else if(mapMode==2) expressionMatrix1 <- expressionMatrix1[,which(colnames(expressionMatrix1) %in% colnames(expressionMatrix2))]
	invisible(expressionMatrix1)
}
