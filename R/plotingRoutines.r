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
plotMapComparison <- function(cross1, cross2, chr=NULL){
	nchrom1 <- length(cross1$geno)
	chrom1 <- as.numeric(names(cross1$geno))
	nchromY <- length(cross2$geno)
	chrom2 <- as.numeric(names(cross2$geno))
	nchrom <- which(chrom1 %in% chrom2)
	cat("Object cross1 contains chromosomes:",paste(chrom1,sep=", "),"and cross2 contains chromosomes:",paste(chrom2,sep=", "),".\n")
	nchromX <- 1:nchrom1
	if(!is.null(chr)){
		nchromX <- NULL
		while(length(chr)>=1){
			if(!(chr[1] %in% nchrom)){
				warning("There are only chromosomes: ",paste(nchrom,sep=", ")," chromosome: ",chr[1]," not found.\n")
			}else{
				nchromX <- c(nchromX,chr[1])
			}
			chr <- chr[-1]
		}
	}
	plotChromosomeMap.internal(cross1,cross2,nchromX,nchromY)

}

compareGeneLocation.internal <- function(cross1, cross2){
	genes1 <- compareGeneLocationSub.internal(cross1)
	print(genes1[1,])
	genes2 <- compareGeneLocationSub.internal(cross2)
	print(genes2[1,])
	genes1 <- mapMarkers.internal(genes1,genes2,mapMode=1)
	genes2 <- mapMarkers.internal(genes2,genes1,mapMode=1)
	result <- matrix(0,nrow(genes1),6)
	rownames(result) <- rownames(genes1)
	for(i in rownames(result)){
		result[i,c(1,2,3)] <- genes1[i,c(1,2,3)]
		result[i,c(4,5,6)] <- genes2[i,c(1,2,3)]
	}
	invisible(result)
}

compareGeneLocationSub.internal <- function(cross){
	genes <- matrix(1,sum(nmar(cross)),3)
	rownames(genes) <- markernames(cross)
	genes[,1] <- rep(1:length(nmar(cross)),nmar(cross))
	genes <- genePosition.internal(genes,cross)
	invisible(genes)
}

genePosition.internal <- function(genes,cross){
	for(i in rownames(genes)){
		chr <- genes[i,1]
		if(chr>1){ 
			prev <- sum(chrlen(cross)[1:(chr-1)]) + (0.13* max(chrlen(cross)) * (chr-1))
		}else{ 
			prev <- 0
		}
		genes[i,2] <- pull.map(cross)[[chr]][i]
		genes[i,3] <- genes[i,2] + prev
	}
	invisible(genes)
}

plotChromosomeMap.internal <- function(cross1,cross2,nchromX,nchromY){
	genes <- compareGeneLocation.internal(cross1,cross2)[which(geneLocationMatrix[,1] %in% nchromX),]
	plot(x=genes[1,3], y=genes[1,6], xlim=c(min(genes[,3]),max(genes[,3])), ylim=c(min(genes[,6]),max(genes[,6])), col="red", xlab="Cross 1", ylab="Cross 2", main="Comparison of genetic maps")
	if(max(nchromY) < max(nchromX)){
		cl <- topo.colors(max(nchromX))
	}else{
		cl <- c(topo.colors(max(nchromX)),rep("red",(max(nchromY) - max(nchromX))))
	}
	for(i in rownames(genes)){
		mark <- 0 + genes[i,1]
		color <- cl[genes[i,4]]
		points(x=genes[i,3], y=genes[i,6],col=color,pch=mark)
	}
	for(x in nchromX){
		if(x > 1) abline(v=sum(chrlen(cross1)[1:(x-1)]) + (0.13* max(chrlen(cross1)) * (x-1)),lty=2)
	}
	for(y in 1:nchromY){
		if(y > 1) abline(h=sum(chrlen(cross2)[1:(y-1)]) + (0.13* max(chrlen(cross2)) * (y-1)),lty=2)
	}
}