############################################################################################################
#
# toGenotypes.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified April 2011
# last modified in version: 0.4.3
# in current version: active, not in main workflow
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
# Contains: toGenotypes
#           convertToGenotypes.internal, splitRow.internal, filterGenotypes.internal, filterRow.internal
#			sortMap.internal, intoRil, removeChromosomes.internal
#
############################################################################################################

############################################################################################################
#cleanMap - removing markers that cause reco map to be too long
# 
# cross - R/qtl cross type object
# difPercentage - If removing marker will make map shorter by this percentage of its length, then it will be dropped.
# verbose - Be verbose
# debugMode - 1: Print our checks, 2: print additional time information 
#
############################################################################################################
cleanMap <- function(cross, difPercentage, verbose=FALSE, debugMode=0){
	if(verbose && debugMode==1) cat("cleanMap starting withour errors in checkpoint.\n")
	s <- proc.time()
	for(i in 1:length(cross$geno)){
		begMarkers <- length(cross$geno[[i]]$map)
		begLength <- max(cross$geno[[i]]$map)
		for(j in names(cross$geno[[i]]$map)){
			if(max(cross$geno[[i]]$map)>150){
				cur_max <- max(cross$geno[[i]]$map)
				cross2 <- drop.markers(cross,j)
				newmap <- est.map(cross2,offset=0)
				cross2 <- replace.map(cross2, newmap)
				new_max <- max(cross2$geno[[i]]$map)
				dif <- cur_max-new_max
				if(dif > (difPercentage/100 * cur_max)){
					if(verbose) cat("------Removed marker:",j,"to make chromosome",i,"map smaller from",cur_max,"to",new_max,"\n")
					cross <- cross2
				}
			}
		}
		removed <- begMarkers-length(cross$geno[[i]]$map)
		if(removed>0)cat("Removed",removed,"out of",begMarkers,"markers on chromosome",i," which led to shortening map from ",begLength,"to",max(cross$geno[[i]]$map),"(",100*(begLength-max(cross$geno[[i]]$map))/begLength,"%)\n")
	}
	e <- proc.time()
	if(verbose && debugMode==2)cat("Map cleaning done in:",(e-s)[3],"seconds.\n")
	invisible(cross)
}

############################################################################################################
#print.ril - overwritting print function for objects of class "ril"
# 
# x - object of class ril
# ... - passed to cats
#
############################################################################################################
print.ril <- function(x,...){
	cat("*************************************************************************************\n")
	cat("This is object of class ril, too complex to print it, so we provide you with summary.\n")
	if(!(is.null(x$rils))){
		if(!(is.null(x$rils$phenotypes))){
			cat("Object contains phenotypic data for",ncol(x$rils$phenotypes),"children individuals covering",nrow(x$rils$phenotypes),"probes.\n",...)
		}else{
			stop("There is no phenotypic data for children in this object. This is not acceptable in real ril object.\n",...)
		}
		if(!(is.null(x$rils$genotypes$read))){
			cat("Object contains genotypic data for",ncol(x$rils$phenotypes),"children individuals covering",nrow(x$rils$phenotypes),"probes.\n",...)
		}else{
			cat("There is no genotypic data for children in this object.\n",...)
		}
		if(!(is.null(x$rils$map))){
			cat("Object contains physical map covering",nrow(x$rils$map),"markers from",length(table(x$rils$map[,1])),"chromosomes.\n",...)
		}else{
			cat("There is no physical genetic map in this object.\n")
		}
	}else{
		cat("WARNING: There is no phenotypic data for children. This is not acceptable in real ril object.\n",...)
	}
	
	if(!(is.null(x$parental))){
		if(!(is.null(x$parental$phenotypes))){
			cat("Object contains phenotypic data for",ncol(x$parental$phenotypes),"parental individuals covering",nrow(x$parental$phenotypes),"probes.\n",...)
		}else{
			stop("There is no phenotypic data for parents in this object. This is not acceptable in real ril object.\n",...)
		}
		if(!(is.null(x$parental$RP))){
			cat("Object contains RP analysis results.\n",...)
		}else{
			cat("There is no RP analysis result in this object.\n",...)
		}
		if(!(is.null(x$parental$groups))){
			cat("Parental groups are as following:",x$parental$groups,"\n",...)
		}else{
			cat("There is no information about parental groups in this object.\n",...)
		}
	}else{
		cat("WARNING: There is no phenotypic data for parents. This is not acceptable in real ril object.\n",...)
	}
	cat("*************************************************************************************\n")
}

############################################################################################################
#intoRil: putting data into ril object
# 
# ril - object of class ril, data should be put into, by default - new object will be created
# parental - matrix of parental data to be put into ril object
# parentalRows - rows to be selected from parental, by default - all
# parentalCols - cols to be selected from parental, by default - all
# children -  matrix of parental data to be put into ril object
# childrenRows - rows to be selected from children, by default - all
# childrenCols - cols to be selected from children, by default - all
#
############################################################################################################
intoRil <- function(ril=NULL, dataObject, dataType=c("parental","children"), selectedRows=NULL, selectedCols=NULL){
	if(!(is.null(dataObject))&&!is.null(dim(dataObject))){
		columns <- NULL
		rows <- NULL
		if(is.null(selectedRows)) selectedRows <- c(1:nrow(dataObject))
		if(is.null(selectedCols)) selectedCols <- c(1:ncol(dataObject))
		### checks
		inRangeCheck.internal(min(selectedRows),"selectedRows",-nrow(dataObject),nrow(dataObject))
		inRangeCheck.internal(max(selectedRows),"selectedRows",-nrow(dataObject),nrow(dataObject))
		inRangeCheck.internal(min(selectedCols),"selectedCols",-ncol(dataObject),ncol(dataObject))
		inRangeCheck.internal(max(selectedCols),"selectedCols",-ncol(dataObject),ncol(dataObject))
		if(length(dataType)!=1) stop("Only one data type could be selected at once.\n")
		inListCheck.internal(dataType,"dataType",c("parental","children"))
		
		dataObject <- dataObject[selectedRows,selectedCols]
		
		for(column in 1:ncol(dataObject)){
			if(!(numericCheck.internal(dataObject[,column],allow.na=TRUE))){
				columns <- c(columns,column)
			}
		}
		
		if(!(is.null(columns))){
			cat("Following  columns are not numeric and cannot be converted into numeric:",columns," so will be removed.\n")
			dataObject <- dataObject[,-columns]
		}
		
		if(is.null(dim(dataObject))) stop("Not enough data to continue.\n")
		
		for(row_ in 1:nrow(dataObject)){
			if(!(numericCheck.internal(dataObject[row_,],allow.na=TRUE))){
				rows <- c(rows,row_)
			}
		}
		
		if(!(is.null(rows))){
			cat("Following  rows are not numeric and cannot be converted into numeric:",rows," so will be removed.\n")
			dataObject <- dataObject[-rows,]
		}
		
		if(is.null(dim(dataObject))) stop("Not enough data to continue.\n")
		
		cur<- matrix(as.numeric(as.matrix(dataObject)),nrow(dataObject),ncol(dataObject))
		
		if(!is.null(rownames(dataObject))){
			rownames(cur) <- rownames(dataObject)
		}else{
			rownames(cur) <- 1:nrow(cur)
		}
		
		if(!is.null(colnames(dataObject))){
			colnames(cur) <- colnames(dataObject)
		}else{
			colnames(cur) <- 1:ncol(cur)
		}
		
		if(dataType=="parental"){
			ril$parental$phenotypes <- cur
			ril$parameters$intoRil$parental <- list("dataObject", dataType, selectedRows, selectedCols)
			names(ril$parameters$intoRil$parental) <- c("dataObject", "dataType", "selectedRows", "selectedCols")
		}else if(dataType=="children"){
			ril$rils$phenotypes <- cur
			ril$parameters$intoRil$children <- list("dataObject", dataType, selectedRows, selectedCols)
			names(ril$parameters$intoRil$children) <- c("dataObject", "dataType", "selectedRows", "selectedCols")
		}
	}
	if(is.null(ril)) stop("No data provided!\n")
	class(ril) <- "ril"
	invisible(ril)
}

############################################################################################################
#removeChromosomes.internal: removing from cross chromosomes that have too little markers
# 
# cross - object of R/qtl cross type
# minChrLength - minimal number of markers chromosome have to contaion not to be removed)
#
############################################################################################################
removeChromosomes.internal <- function(cross, minChrLength){
	j <- length(cross$geno)
	for(i in length(cross$geno):1){
		if(i<=j){
			if(length(cross$geno[[i]]$map)<minChrLength){
				cat("removing markers:",names(cross$geno[[i]]$map),"chr",i,"\n")
				cross$rmv <- cbind(cross$rmv,cross$geno[[i]]$data)
				cross <- drop.markers(cross, names(cross$geno[[i]]$map))
				names(cross$geno) <- 1:length(cross$geno)
				j <- length(cross$geno)
			}
		}
	}
	invisible(cross)
}
