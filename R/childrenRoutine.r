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
# 				mapMarkers, correctChildrenExpression, selectChildrenExpression, correctExpression
#
#################################################################################


#readChildrenExpression
# reading data in, corrects it and uses output from parental function
# to select genes to further analysis

childrenRoutine <- function(childrenFile="Gene_quant.txt",genotypeFile="Genotypes.txt",correction=TRUE,expressionParental,verbose=FALSE,debugMode=0){
	s<-proc.time()
	if(verbose && debugMode==1) cat("readChildrenExpression starting.\n")
	setwd("D:/data/parental")
	library(pheno2geno)
	library(qtl)
	library(iqtl)
	library(RankProd)
	
	genotypeMatrix <- readChildrenGenotypes(genotypeFile)
	expressionChildren <- readExpression(childrenFile,verbose,debugMode)
	expressionChildren <- mapMarkers(expressionChildren,genotypeMatrix)
	genotypeMatrix <- mapMarkers(genotypeMatrix,expressionChildren)
	
	if(correction) expressionChildren <- correctChildrenExpression(expressionChildren,genotypeMatrix,verbose,debugMode)
	expressionChildren <- selectChildrenExpression(expressionChildren,expressionParental,verbose,debugMode)
	
	e<-proc.time()
	if(verbose) cat("readChildrenExpression done in",(e-s)[3],"seconds.\n")
	invisible(expressionChildren)
}

mapMarkers <- function(expressionMatrix1, expressionMatrix2){
	expressionMatrix1 <- expressionMatrix1[,which(colnames(expressionMatrix1) %in% colnames(expressionMatrix2))]
	invisible(expressionMatrix1)
}

correctChildrenExpression <- function(expressionChildren,genotypeMatrix,verbose=FALSE,debugMode=0){
	s2<-proc.time()
	correction <- correctExpression(expressionChildren,genotypeMatrix,verbose,debugMode)
	expressionChildren <- expressionChildren - t(correction)
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
	cross <- genotypesToCross(genotypeMatrix,childrenExpression)
	batchlist <- batcheffectcheck(cross,2,0)
	corrected <- batcheffectcorrect(cross,batchlist,0)
	e<-proc.time()
	if(verbose && debugMode==2)cat("correctExpression done in",(e-s)[3],"seconds.\n")
	invisible(corrected)
}
