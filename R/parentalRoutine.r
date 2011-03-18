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

parentalRoutine <- function(parentalFile="Gene_parental.txt",groupLabels=c(0,0,1,1),treshold=0.01,verbose=FALSE,debugMode=0,...){
	s<-proc.time()
	if(verbose && debugMode==1) cat("readParentalExpression starting.\n")
	setwd("D:/data/parental")
	invisible(library(pheno2geno))
	invisible(library(qtl))
	invisible(library(iqtl))
	invisible(library(RankProd))
	
	expressionParental <- readParentalExpression(parentalFile,verbose,debugMode)
	
	rankParental <- rankParentalExpression(expressionParental,groupLabels,verbose,debugMode,...)
	
	output <- filterParentalExpression(expressionParental,rankParental,groupLabels,treshold,verbose,debugMode)
	
	e<-proc.time()
	if(verbose) cat("readParentalExpression done in",(e-s)[3],"seconds.\n")
	invisible(output)
}

readParentalExpression <- function(parentalFile="Gene_parental.txt",verbose=FALSE,debugMode=0){
	s1<-proc.time()
	expressionParental <- as.matrix(read.table(parentalFile,sep=""))
	e1<-proc.time()
	if(verbose && debugMode==2)cat("Reading parental file:",parentalFile,"done in:",(e1-s1)[3],"seconds.\n")
	invisible(expressionParental)
}

rankParentalExpression <- function(expressionParental,groupLabels=c(0,0,1,1),verbose=FALSE,debugMode=0,...){
	s2<-proc.time()
	rankParental <- invisible(RP(expressionParental,groupLabels,...))
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Product Rank done in:",(e2-s2)[3],"seconds.\n")
	invisible(rankParental)
}

filterParentalExpression <- function(expressionParental,rankParental,groupLabels,treshold=0.01,verbose=FALSE,debugMode=0){
	s2<-proc.time()
	expressionParental <- expressionParental[c(which(rankParental$pval[1]<treshold),which(rankParental$pval[2]<treshold)),]
	output <- matrix(0,nrow(expressionParental),2)
	output[,1] <- apply(expressionParental[,which(groupLabels==0)],1,mean)
	output[,2] <- apply(expressionParental[,which(groupLabels==1)],1,mean)
	rownames(output) <- rownames(expressionParental)
	colnames(output) <- c("Parental_group_0","Parental_group_1")
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Filtering data with treshold:",treshold,"done in:",(e2-s2)[3],"seconds.\n")
	invisible(output)
}