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
#plotParental: plot red points for expression values for parent of type 0, blue for parent 1 and green lines
# for means of rows
#
# ril - Ril type object, must contain parental phenotypic data.
# markers - markers to be printed numbers or names 
# groupLabels - Specify which column of parental data belongs to group 0 and which to group 1
#
############################################################################################################
plotParentalExpression <- function(ril, markers=1:100, groupLabels=c(0,0,1,1)){
	if(nrow(ril$parental$phenotypes)<length(markers)){
		cat("WARNING: There are only",nrow(ril$parental$phenotypes),"markers in ril$parental$phenotypes, function will plot them all.\n")
		markers <- 1:nrow(ril$parental$phenotypes)
	}
	for(i in markers) if(!(i %in% 1:nrow(ril$parental$phenotypes))){
		stop("ERROR: There is no marker number: ",i,", stopping \n")
	}
	parental <- ril$parental$phenotypes[markers,]
	plot(x=markers[1], y=parental[1,1], xlim=c(min(markers),max(markers)), ylim=c(min(parental),max(parental)), col="red",
	xlab="Marker", ylab="Expression value", main="Parental gene expression data")
	for(i in 1:nrow(parental)){
		for(j in which(groupLabels==0)){
			if(i%%50==0) cat(i,j,parental[i,j],"\n")
			points(x=markers[i],y=parental[i,j],col="red")
		}
		for(k in which(groupLabels==1)){
			if(i%%50==0) cat(i,k,parental[i,k],"\n")
			points(x=markers[i],y=parental[i,k],col="blue")
		}
	}
	points(x=markers,y=apply(parental,1,mean),col="green", pch=95, cex=3)
}

############################################################################################################
#plotChildrenExpression: boxplot of data for selected markers + points of parental mean for each marker
# 
# ril - Ril type object, must contain parental phenotypic data.
# markers - markers to be printed numbers or names 
#
############################################################################################################
plotChildrenExpression <- function(ril, markers=1:100){
	if(nrow(ril$rils$phenotypes)<length(markers)){
		cat("WARNING: There are only",nrow(ril$rils$phenotypes),"markers in ril$rils$phenotypes, function will plot them all.\n")
		markers <- 1:nrow(ril$rils$phenotypes)
	}
	if(nrow(ril$parental$phenotypes)<length(markers)){
		cat("WARNING: There are only",nrow(ril$parental$phenotypes),"markers in ril$parental$phenotypes, function will plot them all.\n")
		markers <- 1:nrow(ril$parental$phenotypes)
	}
	for(i in markers) if(!(i %in% 1:nrow(ril$rils$phenotypes))){
		stop("ERROR: There is no marker number: ",i,", stopping \n")
	}
	for(i in markers) if(!(i %in% 1:nrow(ril$parental$phenotypes))){
		stop("ERROR: There is no marker number: ",i,", stopping \n")
	}
	children <- ril$rils$phenotypes[markers,]
	parental <- ril$parental$phenotypes[which(rownames(ril$parental$phenotypes) %in% rownames(children)),]
	rownames(children) <- markers
	boxplot(t(children), ylim=c(min(children), max(children)),	xlab="Marker", ylab="Expression value", main="Children gene expression data")
	points(apply(parental,1,mean),col="green", pch=95, cex=3)
	points(apply(parental,1,max),col="red", pch=24, cex=1)
	points(apply(parental,1,min),col="blue", pch=25, cex=1)
}

############################################################################################################
#plotChildrenExpression: boxplot of data for selected markers + points of parental mean for each marker
# 
# ril - Ril type object, must contain parental phenotypic data.
# markers - markers to be printed numbers or names 
#
############################################################################################################
plotMapComparison <- function(cross){
	referenceLocs <- getXlocs.internal(cross)
	predictedLocs <- getYLocs.internal(cross)
	plot(x=referenceLocs[[1]], y=predictedLocs[[1]], xlim=c(min(x=referenceLocs[[1]]),max(x=referenceLocs[[1]]), ylim=c(min(predictedLocs[[1]]),max(predictedLocs[[1]])),
		col="red", xlab="Reference map", ylab="Predicted map", main="Comparison of genetic maps")
	for(x in 2:length(referenceLocs[[2]])){
		abline(v=sum(referenceLocs[[2]][1:(x-1)]) + (0.13* max(referenceLocs[[2]]) * (x-1)),lty=2)
	}
	for(y in 2:length(predictedLocs[[2]])){
		abline(h=sum(predictedLocs[[2]][1:y-1]) + (0.13* max(predictedLocs[[2]]) * (y-1)),lty=2)
	}
}

getChromosome.internal <- function(cross){
	invisible(rep(1:length(nmar(cross)),nmar(cross)))
}

getYLocs.internal <- function(cross){
	locs <- lapply(cross$geno,function(x){as.numeric(x$map)})
	chrlength <- lapply(locs,max)
	summaryLengths <- rep(0,length(chrlength))
	for(x in 2:length(chrlength)){
		summaryLengths[x] <- max(chrlength[[x-1]]) + summaryLengths[x-1] + 1000
	}
	chrids <- getChromosome(cross)
	result <- as.matrix(summaryLengths[chrids]+as.numeric(unlist(locs)))
	rownames(result) <- markernames(cross)
	invisible(list(result,chrlength))
}

getXlocs.internal (cross,markers){
	cross$maps$physical 
}