#################################################################################
#
# dataRoutines.R
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
# Contains: readParentalExpression, readChildrenExpression, correctExpression
#           
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

readParentalExpression <- function(parentalFile="Gene_parental.txt",groupLabels=c(0,0,1,1),treshold=0.01,verbose=FALSE,debugMode=0,...){
	s<-proc.time()
	if(verbose && debugMode==1) cat("recombinationCount starting.\n")
	setwd("D:/data/parental")
	library(pheno2geno)
	library(qtl)
	library(iqtl)
	library(RankProd)
	
	s1<-proc.time()
	expressionParental <- as.matrix(read.table(parentalFile,sep=""))
	e1<-proc.time()
	if(verbose && debugMode==2)cat("Reading parental file:",parentalFile,"done in:",(e1-s1)[3],"seconds.\n")
	
	s2<-proc.time()
	rankParental <- invisible(RP(expressionParental,groupLabels,...))
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Product Rank done in:",(e2-s2)[3],"seconds.\n")
	
	expressionParental <- expressionParental[c(which(rankParental$pval[1]<treshold),which(rankParental$pval[2]<treshold)),]
	output <- matrix(0,nrow(expressionParental),2)
	output[,1] <- apply(expressionParental[,which(groupLabels==0)],1,mean)
	output[,2] <- apply(expressionParental[,which(groupLabels==1)],1,mean)
	rownames(output) <- rownames(expressionParental)
	colnames(output) <- c("Parental_group_0","Parental_group_1")
	
	e<-proc.time()
	if(verbose) cat("recombinationCount done in",(e-s)[3],"seconds.\n")
	invisible(output)
}


readChildrenExpression <- function(filename,parentalOutput, correction=TRUE){
	#reading data in, corrects it and uses output from parental function
	# to select genes to further analysis
}

readChildrenGenotypes <- function(filename="Genotypes.txt"){
	res1 <- as.matrix(read.table(filename,sep="\t",header=T))
	res1 <- t(res1)
	res2 <- res1[-1,]
	rownames(res2)<-rownames(res1[-1,], do.NULL = FALSE)
	colnames(res2)<-res1[1,]
	invisible(res2)
}


correctExpression <- function(){
	#using Danislav's batchcorrection probably, must think if it is wise to produce uber-massive
	#cross file...
}
