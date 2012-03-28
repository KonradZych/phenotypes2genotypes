#################################################################################
#
# find.diff.expressed.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified November 2011
# last modified in version: 0.9.1
# in current version: active, in main workflow
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
# Contains: find.diff.expressed
#				fake.population, plotRPpval
#
#################################################################################

############################################################################################################
#									*** find.diff.expressed ***
#
# DESCRIPTION:
#	Using Rank Product or student t-test analysis to select differentially expressed genes.
# 
# PARAMETERS:
# 	population - Object of class population , must contain founders phenotypic data.
# 	verbose - Be verbose
# 	debugMode - 1: Print our checks, 2: print additional time information
# 	... - parameters send to RP function
# 
# OUTPUT:
#	object of class population containing object of class RP in $founders$RP
#
############################################################################################################
find.diff.expressed <- function(population,use=c("ttest","rankprod"),verbose=FALSE,debugMode=0,...){
  #checks
  if(missing(population)) stop("provide population object\n")
  check.population(population)
  use <- checkParameters.internal(use,c("ttest","rankprod"),"use")
  if(verbose && debugMode==1) cat("find.diff.expressed starting withour errors in checkpoints.\n")
	
	s<-proc.time()
	if(use=="rankprod"){
		tryCatch(require(RankProd), error = stop("Install RankProd package to use Rank Product analysis!\n"))
		rankProdRes <- RP(population$founders$phenotypes,population$founders$groups,gene.names=rownames(population$founders$phenotypes),...)
		population$founders$RP <- rankProdRes
	}else{
		population$founders$RP$pval<- t(rbind(apply(population$founders$phenotypes,1,findUsingTTest.internal,population$founders$groups)))
	}
	e<-proc.time()
	if(verbose && debugMode==2)cat("Differentially expressed genes found in:",(e-s)[3],"seconds.\n")
	invisible(population)
}

############################################################################################################
#									*** findUsingTTest.internal ***
#
# DESCRIPTION:
#	subfunction of find.diff.expressed using t-test to assess whether gene is differentially expressed
# 
# PARAMETERS:
# 	phenoRow - single row of founders phenotype data
# 	groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
#
# OUTPUT:
#	two p-values - for gene being up- and downregulated
#
############################################################################################################
findUsingTTest.internal <- function(phenoRow,groupLabels){
  a <- which(groupLabels==0)
  b <- which(groupLabels==1)
  if(mean(phenoRow[a]) < mean(phenoRow[b])){
    what <- "less"
    return(c(0,t.test(phenoRow[a],phenoRow[b],alt=what)$p.value))
  }else{
    what <- "gre"
    return(c(t.test(phenoRow[a],phenoRow[b],alt=what)$p.value,0))
  }
}

############################################################################################################
#									*** showRPpval ***
#
# DESCRIPTION:
#	showing pvals of RP for selected markers
# 
# PARAMETERS:
# 	population - Object of class population , must contain founders phenotypic data.
# 	markers - markers (specified by number) to be shown
# 
# OUTPUT:
#	none
#
############################################################################################################
showRPpval <- function(population,markers=1:10){
	#checks
  if(missing(population)) stop("provide population object\n")
  check.population(population)
  if(is.null(population$founders$RP$pval)) stop("Population object does not contain results of RP analysis run find.diff.expressed first.\n")
  inRangeCheck.internal(markers,"markers",1,nrow(population$founders$phenotypes))
  
	toPrint <- matrix(0,length(markers),2)
	toPrint[,1] <- population$founders$RP$pval[markers,1]
	toPrint[,2] <- population$founders$RP$pval[markers,2]
	rownames(toPrint) <- rownames(population$founders$phenotypes)[markers]
	colnames(toPrint) <- c("up","down")
	print(toPrint)
}

############################################################################################################
#									*** plotRPpval ***
#
# DESCRIPTION:
#	ploting pvals of RP for selected markers
# 
# PARAMETERS:
# 	population - Object of class population , must contain founders phenotypic data.
# 	markers - markers (specified by number) to be shown
#	treshold - treshold value, on which line is plotted (by default - 0.01)
# 
# OUTPUT:
#	none
#
############################################################################################################
plotRPpval <- function(population,thresholdRange=c(0.01,0.1,0.01)){
	#checks
  if(missing(population)) stop("provide population object\n")
  check.population(population)
  if(is.null(population$founders$RP$pval)) stop("Population object does not contain results of RP analysis run find.diff.expressed first.\n")
  thrRange <- seq(thresholdRange[1],thresholdRange[2],thresholdRange[3])
  n.upSelected <- NULL
  n.downSelected <- NULL
  for(threshold in thrRange){
    upNotNull <- population$founders$RP$pval[which(population$founders$RP$pval[,1] > 0),1]
    downNotNull <- population$founders$RP$pval[which(population$founders$RP$pval[,2] > 0),2]
    n.upSelected <- c(n.upSelected,length(which(upNotNull < threshold)))
    n.downSelected <- c(n.downSelected,length(which(downNotNull < threshold)))
  }
	plot(thrRange,n.upSelected ,main="RP analysis p-values",xlab="p-value",ylab="# markers selected",xlim=c(thresholdRange[1],thresholdRange[2]),ylim=c(min(min(n.upSelected),min(n.downSelected )),max(max(n.upSelected),max(n.downSelected))),type="o")
  points(thrRange,n.downSelected,col="red",type="o")
  legend(x="topleft",legend=c("up regulated","down regulated"),col=c("black","red"),cex=0.8,pch=21,lwd=2,bg="white")
}
