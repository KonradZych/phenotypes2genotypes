##################################################################################################
#
# plottingRoutines.R
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
# Contains: 
#           
#
##################################################################################################


plotParental <- function(parentalExpression,groupLabels=c(0,1),plotMode="rank",rankParental=NULL,treshold=0.05,verbose=FALSE,debugMode=0){
	if(plotMode=="rank"){
		colPal <- matrix("gray",nrow(parentalExpression),ncol(parentalExpression))
		colPal[which(rankParental$pval[1]<treshold),] <- "red"
		colPal[which(rankParental$pval[2]<treshold),] <- "black"
	}else{
		colPal <- matrix("black",nrow(parentalExpression),ncol(parentalExpression))
		colPal[which(groupLabels==1),] <- "red"
	}
	plot(parentalExpression,col=colPal)
}