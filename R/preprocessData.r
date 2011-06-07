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
# population - Population type object, must contain founders phenotypic data.
# groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information
#
############################################################################################################
preprocessData <- function(population,groupLabels=c(0,0,1,1),verbose=FALSE,debugMode=0,...){
	s2<-proc.time()
	rankProdRes <- cat(RP(population$founders$phenotypes,groupLabels,...),file=NULL)
	population$founders$RP <- rankProdRes
	population$founders$groups <- groupLabels
	population$parameters$preprocessData <- list("object of population class", groupLabels,verbose,debugMode)
	names(population$parameters$preprocessData) <- c("population", "groupLabels", "verbose", "debugMode")
	class(population) <- "population"
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Data preprocessing done in:",(e2-s2)[3],"seconds.\n")
	invisible(population)
}
