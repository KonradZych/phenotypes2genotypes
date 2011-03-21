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

childrenRoutine <- function(){
	s<-proc.time()
	cat("childrenRoutine starting.\n")
	setwd("D:/data/parental")
	invisible(require(iqtl))
	invisible(require(pheno2geno))
	invisible(require(RankProd))
	
	genotypeMatrix <- readChildrenGenotypes("Genotypes.txt",verbose=TRUE,debugMode=2)
	expressionChildren <- readExpression("Gene_quant.txt",verbose=TRUE,debugMode=2)
	expressionChildren <- mapMarkers(expressionChildren,genotypeMatrix)
	genotypeMatrix <- mapMarkers(genotypeMatrix,expressionChildren)
	
	expressionChildren <- correctChildrenExpression(expressionChildren,genotypeMatrix,verbose=TRUE,debugMode=2)
	expressionChildrenSel <- selectChildrenExpression(expressionChildren,expressionParental,verbose=TRUE,debugMode=2)
	
	e<-proc.time()
	cat("readChildrenExpression done in",(e-s)[3],"seconds.\n")
	invisible(list(expressionChildren,expressionChildrenSel[[1]],expressionChildrenSel[[2]]))
}

mapMarkers <- function(expressionMatrix1, expressionMatrix2, mapMode=2){
	if(mapMode==1) expressionMatrix1 <- expressionMatrix1[which(rownames(expressionMatrix1) %in% rownames(expressionMatrix2)),]
	else if(mapMode==2) expressionMatrix1 <- expressionMatrix1[,which(colnames(expressionMatrix1) %in% colnames(expressionMatrix2))]
	invisible(expressionMatrix1)
}

correctChildrenExpression <- function(expressionChildren,genotypeMatrix,verbose=FALSE,debugMode=0){
	s2<-proc.time()
	if(verbose && debugMode==1) cat("correctChildrenExpression starting.\n")
	correction <- correctExpression(expressionChildren,genotypeMatrix,verbose,debugMode)
	expressionChildren <- expressionChildren - (correction)
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Correcting expression data done in:",(e2-s2)[3],"seconds.\n")
	invisible(expressionChildren)
}

selectChildrenExpression <- function(expressionChildren,expressionParental,verbose=FALSE,debugMode=0){
	s2<-proc.time()
	if(verbose && debugMode==1) cat("selectChildrenExpression starting.\n")
	expressionChildrenUp <- expressionChildren[which(rownames(expressionChildren) %in% expressionParental[[2]]),]
	expressionChildrenDown <- expressionChildren[which(rownames(expressionChildren) %in% expressionParental[[3]]),]
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Selecting expression data done in:",(e2-s2)[3],"seconds.\n")
	invisible(list(expressionChildrenUp,expressionChildrenDown))
}

readChildrenGenotypes <- function(filename,verbose=FALSE,debugMode=0){
	s1<-proc.time()
	res1 <- as.matrix(read.table(filename,sep="\t",header=T))
	res1 <- t(res1)
	res2 <- matrix(as.character("-"),nrow(res1)-1,ncol(res1))
	res2[which(res1[-1,]=="A")] <- 0
	res2[which(res1[-1,]=="a")] <- 0
	res2[which(res1[-1,]=="B")] <- 1
	res2[which(res1[-1,]=="b")] <- 1
	rownames(res2)<-rownames(res1[-1,], do.NULL = FALSE)
	colnames(res2)<-res1[1,]
	e1<-proc.time()
	if(verbose && debugMode==2)cat("Reading genotypes file:",filename,"done in:",(e1-s1)[3],"seconds.\n")
	invisible(res2)
}

correctExpression <- function(expressionMatrix,genotypeMatrix,verbose=FALSE,debugMode=0){
	s<-proc.time()
	if(verbose && debugMode==1) cat("correctExpression starting.\n")
	if(file.exists("tmpbatch.out")){
		if(verbose) cat("File tmpbatch.out already exists, reading it.\n")
		s1 <- proc.time()
		corrected <- read.table("tmpbatch.out",sep="\t")
		corrected <- corrected[1:nrow(expressionMatrix),1:ncol(expressionMatrix)]
		e1 <- proc.time()
		if(verbose && debugMode==2)cat("Reading tmpbatch.out done in:",(e1-s1)[3],"seconds.\n")
	}else{
		if(file.exists("realgenotypescross.csv")){
			if(verbose) cat("File realgenotypescross.csv already exists, reading it.\n")
			s2 <- proc.time()
			cross <- invisible(read.cross("csvr",file="realgenotypescross.csv", genotypes=c(0,1)))
			e2 <- proc.time()
			if(verbose && debugMode==2)cat("Reading realgenotypescross.csv done in:",(e2-s2)[3],"seconds.\n")
		}else{
			cross <- genotypesToCross(genotypeMatrix,expressionMatrix,outputFile="realgenotypescross.csv",verbose=verbose,debugMode=debugMode)
		}
		batchlist <- batcheffectcheck(cross,2,0)
		corrected <- t(batcheffectcorrect(cross,batchlist,0))
	}
	e<-proc.time()
	if(verbose && debugMode==2)cat("correctExpression done in",(e-s)[3],"seconds.\n")
	invisible(corrected)
}
