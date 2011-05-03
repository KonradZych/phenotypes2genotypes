#################################################################################
#
# preprocessData.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified May 2011
# last modified in version: 0.5.1
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
#toGenotypes: Using Rank Product analysis to select differentially expressed genes.
# 
# ril - Ril type object, must contain parental phenotypic data.
# groupLabels - Specify which column of parental data belongs to group 0 and which to group 1.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
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
	ril$parental$groups <- groupLabels
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Data preprocessing done in:",(e2-s2)[3],"seconds.\n")
	ril$parameters$preprocessData <- list("object of ril class", groupLabels,verbose,debugMode)
	names(ril$parameters$preprocessData) <- c("ril", "groupLabels", "verbose", "debugMode")
	class(ril) <- "ril"
	invisible(ril)
}

