#################################################################################
#
# childrenRoutine.R
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
# Contains: childrenRoutine, readChildrenGenotypes, 
# 				readChildrenExpression, mapMarkers, correctChildrenExpression, selectChildrenExpression, correctExpression
#
#################################################################################


#readChildrenExpression
# reading data in, corrects it and uses output from parental function
# to select genes to further analysis

childrenRoutine <- function(childrenFile="Gene_quant.txt",genotypeFile="Genotypes.txt",correction=TRUE,expressionParental=NULL,verbose=FALSE,debugMode=0){
	s<-proc.time()
	if(verbose && debugMode==1) cat("readChildrenExpression starting.\n")
	setwd("D:/data/parental")
	library(pheno2geno)
	library(qtl)
	library(iqtl)
	library(RankProd)
	
	genotypeMatrix <- readChildrenGenotypes(genotypeFile)
	childrenExpression <- readChildrenExpression(childrenFile,verbose,debugMode)
	childrenExpression <- mapMarkers(childrenExpression,genotypeMatrix)
	print(dim(childrenExpression))
	genotypeMatrix <- mapMarkers(genotypeMatrix,childrenExpression)
	print(dim(genotypeMatrix))
	
	if(correction) childrenExpression <- correctChildrenExpression(childrenExpression,genotypeMatrix,verbose,debugMode)
	if(expressionParental!=NULL) childrenExpression <- selectChildrenExpression(childrenExpression,expressionParental,verbose,debugMode)
	
	e<-proc.time()
	if(verbose) cat("readChildrenExpression done in",(e-s)[3],"seconds.\n")
	invisible(expressionChildren)
}

readChildrenExpression <- function(childrenFile,verbose=FALSE,debugMode=0){
	s1<-proc.time()
	expressionChildren <- as.matrix(read.table(childrenFile,sep=""))
	e1<-proc.time()
	if(verbose && debugMode==2)cat("Reading children file:",childrenFile,"done in:",(e1-s1)[3],"seconds.\n")
	invisible(expressionChildren)
}

mapMarkers <- function(expressionMatrix1, expressionMatrix2){
	expressionMatrix1 <- expressionMatrix1[,which(colnames(expressionMatrix1) %in% colnames(expressionMatrix2))]
	invisible(expressionMatrix1)
}

correctChildrenExpression <- function(childrenExpression,genotypeMatrix,verbose=FALSE,debugMode=0){
	s2<-proc.time()
	expressionChildren <- expressionChildren + t(correctExpression(expressionChildren,genotypeMatrix,verbose,debugMode))
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Correcting expression data done in:",(e2-s2)[3],"seconds.\n")
	invisible(expressionChildren)
}

selectChildrenExpression <- function(expressionChildren,expressionParental,verbose=FALSE,debugMode=0){
	s2<-proc.time()
	expressionChildren <- expressionChildren[which(rownames(expressionParental) %in% rownames(expressionChildren)),]
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Selecting expression data done in:",(e2-s2)[3],"seconds.\n")
	invisible(expressionChildren)
}

readChildrenGenotypes <- function(filename){
	res1 <- as.matrix(read.table(filename,sep="\t",header=T))
	res1 <- t(res1)
	res2 <- matrix(as.character("-"),nrow(res1)-1,ncol(res1))
	res2[which(res1[-1,]=="A")] <- 0
	res2[which(res1[-1,]=="a")] <- 0
	res2[which(res1[-1,]=="B")] <- 1
	res2[which(res1[-1,]=="b")] <- 1
	rownames(res2)<-rownames(res1[-1,], do.NULL = FALSE)
	colnames(res2)<-res1[1,]
	invisible(res2)
}

correctExpression <- function(expressionMatrix,genotypeMatrix,verbose=FALSE,debugMode=0){
	s<-proc.time()
	if(verbose && debugMode==1) cat("correctExpression starting.\n")
	cross <- genotypesToCross(genotypeMatrix,expressionMatrix,verbose=verbose,debugMode=debugMode)
	batchlist <- batcheffectcheck(cross,2,0)
	corrected <- batcheffectcorrect(cross,batchlist,0)
	e<-proc.time()
	if(verbose && debugMode==2)cat("correctExpression done in",(e-s)[3],"seconds.\n")
	invisible(corrected)
}
