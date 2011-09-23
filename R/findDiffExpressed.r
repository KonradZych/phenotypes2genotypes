#################################################################################
#
# findDiffExpressed.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified September 2011
# last modified in version: 0.9.0
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
# Contains: findDiffExpressed
#				showRPpval, plotRPpval
#
#################################################################################

############################################################################################################
#									*** findDiffExpressed ***
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
findDiffExpressed <- function(population,use=c("ttest","rankprod"),verbose=FALSE,debugMode=0,...){
	if(missing(population)) stop("provide population object\n")
	use <- defaultCheck.internal(use,"use",2,"ttest")
	is.population(population)
	s<-proc.time()
	if(use=="rankprod"){
		rankProdRes <- RP(population$founders$phenotypes,population$founders$groups,gene.names=rownames(population$founders$phenotypes),...)
		population$founders$RP <- rankProdRes
	}else{
		population$founders$RP$pval<- t(rbind(apply(population$founders$phenotypes,1,findUsingTTest.internal,population$founders$groups)))
	}
	class(population) <- "population"
	e<-proc.time()
	if(verbose && debugMode==2)cat("Differentially expressed genes found in:",(e-s)[3],"seconds.\n")
	invisible(population)
}

############################################################################################################
#									*** findUsingTTest.internal ***
#
# DESCRIPTION:
#	subfunction of findDiffExpressed using t-test to assess whether gene is differentially expressed
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
ttest_res <- t.test(phenoRow[a],phenoRow[b])
ttest_res <- unlist(ttest_res)
mean_x <- as.numeric(ttest_res[6])
mean_y <- as.numeric(ttest_res[7])
if(mean_x<mean_y){
	return(c(as.numeric(ttest_res[3]),1-(as.numeric(ttest_res[3]))))
}else{
	return(c(1-as.numeric(ttest_res[3]),(as.numeric(ttest_res[3]))))
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
	if(missing(population)) stop("provide population object\n")
	if(min(markers<1)||max(markers)>nrow(population$founders$phenotypes)) stop("wrong range of markers selected\n")
	is.population(population)
	if(is.null(population$founders$RP$pval)) stop("population object does not contain results of RP analysis\n")
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
plotRPpval <- function(population,markers=1:10,treshold=0.01){
	if(missing(population)) stop("provide population object\n")
	if(min(markers<1)||max(markers)>nrow(population$founders$phenotypes)) stop("wrong range of markers selected\n")
	is.population(population)
	if(is.null(population$founders$RP$pval)) stop("population object does not contain results of RP analysis\n")
	plot(population$founders$RP$pval[markers,1],main="RP analysis p-values",xlab="markers",ylab="p-value",ylim=c(0,1))
	points(population$founders$RP$pval[markers,2],col="red")
	abline(h=treshold)
}
