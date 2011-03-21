#################################################################################
#
# parentalRoutine.R
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
# Contains: parentalRoutine 
# 				readParentalExpression, rankParentalExpression, filterParentalExpression
#
#################################################################################

#readParentalExpression
# reading data in, user need to specify groups - sequence of 0s and 1s
# this will be used in product rank function -> selected genes are loaded
# parentalFile
# groupLabels
# treshold
# verbose - standard
# debugMode - standard

parentalRoutine <- function(){
	s<-proc.time()
	cat("readParentalExpression starting.\n")
	setwd("D:/data/parental")
	invisible(require(iqtl))
	invisible(require(pheno2geno))
	invisible(require(RankProd))

	expressionParental <- readExpression("Gene_parental.txt",verbose=TRUE,debugMode=2)

	rankParental <- rankParentalExpression(expressionParental,groupLabels=c(0,0,1,1),verbose=TRUE,debugMode=2)

	parental <- filterParentalExpression(expressionParental,rankParental,treshold=0.01,verbose=TRUE,debugMode=2)

	e<-proc.time()
	cat("readParentalExpression done in",(e-s)[3],"seconds.\n")
	invisible(list(expressionParental,parental[[1]],parental[[2]]))
}

readExpression <- function(expressionFile="Gene_parental.txt",verbose=FALSE,debugMode=0){
	s1<-proc.time()
	expressionMatrix <- as.matrix(read.table(expressionFile,sep=""))
	e1<-proc.time()
	if(verbose && debugMode==2)cat("Reading expression file:",expressionFile,"done in:",(e1-s1)[3],"seconds.\n")
	invisible(expressionMatrix)
}

rankParentalExpression <- function(expressionParental,groupLabels=c(0,0,1,1),verbose=FALSE,debugMode=0,...){
	s2<-proc.time()
	if(file.exists("RP.Rdata")){
		if(verbose) cat("File RP.Rdata already exists, reading it.\n")
		load("RP.Rdata")
	}else{
		rankParental <- invisible(RP(expressionParental,groupLabels,...))
		save(file="RP.Rdata",rankParental)
	}
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Product Rank done in:",(e2-s2)[3],"seconds.\n")
	invisible(rankParental)
}

filterParentalExpression <- function(expressionParental,rankParental,treshold=0.01,verbose=FALSE,debugMode=0){
	s2<-proc.time()
	up <- rownames(expressionParental[which(rankParental$pval[1] < treshold),])
	down <- rownames(expressionParental[which(rankParental$pval[2] < treshold),])
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Filtering data with treshold:",treshold,"done in:",(e2-s2)[3],"seconds.\n")
	invisible(list(up,down))
}