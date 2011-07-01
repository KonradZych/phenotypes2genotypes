#################################################################################
#
# findDiffExpressed.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified July 2011
# last modified in version: 0.8.1
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
# Contains: findDiffExpressed 
#
#################################################################################

############################################################################################################
#									*** findDiffExpressed ***
#
# DESCRIPTION:
#	Using Rank Product analysis to select differentially expressed genes.
# 
# PARAMETERS:
# 	population - Object of class population , must contain founders phenotypic data.
# 	groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
# 	verbose - Be verbose
# 	debugMode - 1: Print our checks, 2: print additional time information
# 	... - parameters send to RP function
# 
# OUTPUT:
#	object of class population containing object of class RP in $founders$RP
#
############################################################################################################
findDiffExpressed <- function(population,verbose=FALSE,debugMode=0,...){
	is.population(population)
	s<-proc.time()
	if(is.null(population$founders$phenotypes)) stop("No founders phenotype data provided\n")
	if(is.null(population$founders$groups)) stop("No information about founders groups data provided\n")
	rankProdRes <- invisible(RP(population$founders$phenotypes,population$founders$groups,...))
	population$founders$RP <- rankProdRes
	class(population) <- "population"
	e<-proc.time()
	if(verbose && debugMode==2)cat("Data preprocessing done in:",(e-s)[3],"seconds.\n")
	invisible(population)
}
