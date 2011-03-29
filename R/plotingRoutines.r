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


############################################################################################################
#plotParental: not that useful, but stays for now
# 
############################################################################################################
plotParental <- function(ril, markers=1:100, groupLabels=c(0,0,1,1), verbose=FALSE, debugMode=0){
	stop("due to stupid R bug, no you shouldn't use it\n")
	parental <- ril$parental$phenotypes[markers,]
	# R BUG - really stupid, xlim should be (x1,x2), but if I put it this way I get error that there is "," that
	# is not allowed...
	#plot(x=markers[1], y=parental[which(groupLabels==0)[1]], xlim=(1,length(markers)), col="red" )
	for(i in 1:length(markers)){
		for(j in 1:ncol(parental[,which(groupLabels==0)])){
			points(x=i,y=parental[j],col="red")
		}
		for(j in 1:ncol(parental[,which(groupLabels==1)])){
			points(x=i,y=parental[j],col="blue")
		}
	}
	points(apply(parental,1,mean),col="green", pch=95, cex=3)
}

############################################################################################################
#plotChildrenExpression: boxplot of data for selected markers + points of parental mean for each marker
# 
# ril - Ril type object, must contain parental phenotypic data.
# markers - markers to be printed numbers or names 
# ... - passed to boxplot
############################################################################################################
plotChildrenExpression <- function(ril, markers=1:100, ...){
	boxplot(t(ril$rils$phenotypes[markers,]))
	parental <- ril$parental$phenotypes[which(rownames(ril$parental$phenotypes) %in% rownames(ril$rils$phenotypes[markers,])),]
	points(apply(parental,1,mean),col="green", pch=95, cex=3)
	points(apply(parental,1,max),col="blue", pch=24, cex=1)
	points(apply(parental,1,min),col="red", pch=25, cex=1)
}