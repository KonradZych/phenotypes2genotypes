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
#print.population - overwritting print function for objects of class "population"
# 
# x - object of class population
# ... - passed to cats
#
############################################################################################################
print.population <- function(x,...){
	cat("*************************************************************************************\n")
	cat("This is object of class population, too complex to print it, so we provide you with summary.\n")
	if(!(is.null(x$offspring))){
		if(!(is.null(x$offspring$phenotypes))){
			cat("Object contains phenotypic data for",ncol(x$offspring$phenotypes),"offspring individuals covering",nrow(x$offspring$phenotypes),"probes.\n",...)
		}else{
			stop("There is no phenotypic data for offspring in this object. This is not acceptable in real ril object.\n",...)
		}
		if(!(is.null(x$offspring$genotypes$read))){
			cat("Object contains genotypic data for",ncol(x$offspring$phenotypes),"offspring individuals covering",nrow(x$offspring$phenotypes),"probes.\n",...)
		}else{
			cat("There is no genotypic data for offspring in this object.\n",...)
		}
		if(!(is.null(x$offspring$map))){
			cat("Object contains physical map covering",nrow(x$offspring$map),"markers from",length(table(x$offspring$map[,1])),"chromosomes.\n",...)
		}else{
			cat("There is no physical genetic map in this object.\n")
		}
	}else{
		cat("WARNING: There is no phenotypic data for offspring. This is not acceptable in real ril object.\n",...)
	}
	
	if(!(is.null(x$founders))){
		if(!(is.null(x$founders$phenotypes))){
			cat("Object contains phenotypic data for",ncol(x$founders$phenotypes),"founders individuals covering",nrow(x$founders$phenotypes),"probes.\n",...)
		}else{
			stop("There is no phenotypic data for parents in this object. This is not acceptable in real ril object.\n",...)
		}
		if(!(is.null(x$founders$RP))){
			cat("Object contains RP analysis results.\n",...)
		}else{
			cat("There is no RP analysis result in this object.\n",...)
		}
		if(!(is.null(x$founders$groups))){
			cat("Parental groups are as following:",x$founders$groups,"\n",...)
		}else{
			cat("There is no information about founders groups in this object.\n",...)
		}
	}else{
		cat("WARNING: There is no phenotypic data for parents. This is not acceptable in real ril object.\n",...)
	}
	cat("*************************************************************************************\n")
}

############################################################################################################
#createPopulation: creating population object
# 
# offspring$phenotypes - matrix containing offspring phenotype data (have to be supported, if not - function
#	quits with error
# founders - matrix containing founders phenotype data (optional)
# offspring$genotypes - matrix containing offspring genotype data (optional)
# maps$genetic - matrix containing genetic map (optional)
# maps$physical - matrix containing physical map (optional)
#
############################################################################################################
createPopulation <- function(offspring_phenotypes=NULL, founders=NULL, offspring_genotypes=NULL, maps_genetic=NULL, maps_physical=NULL,verbose=FALSE,debugMode=0){
	if(verbose && debugMode==1) cat("createPopulation starting.\n")
	s <- proc.time()
	population <- NULL
	if(is.null(offspring_phenotypes)){
		stop("No offspring phenotype data provided!\n")
	}else{
		population <- intoPopulationSub.internal(population, offspring_phenotypes, "offspring$phenotypes", verbose, debugMode)
	}
	if(is.null(founders)){
		warning("No founders phenotype data provided. Strongly recommend to supply this data using intoPopulation.\n")
	}else{
		population <- intoPopulationSub.internal(population, founders, "founders", verbose, debugMode)
	}
	if(is.null(offspring_genotypes)){
		if(verbose)cat("No offspring genotypic data provided. You can supply it later using intoPopulation.\n")
	}else{
		population <- intoPopulationSub.internal(population, offspring_genotypes, "offspring$genotypes", verbose, debugMode)
	}
	if(is.null(maps_genetic)){
		if(verbose)cat("No genotic map provided. You can supply it later using intoPopulation.\n")
	}else{
		population <- intoPopulationSub.internal(population, maps_genetic, "maps$genetic", verbose, debugMode)
	}
	if(is.null(offspring_phenotypes)){
		if(verbose)cat("No physical map provided.  You can supply it later using intoPopulation.\n")
	}else{
		population <- intoPopulationSub.internal(population, maps_physical, "maps$physical", verbose, debugMode)
	}
	if(is.null(population)) stop("No data provided!\n")
	class(population) <- "population"
	e <- proc.time()
	cat("createPopulation done in:",(e-s)[3],"seconds.\n")
	invisible(population)
}

############################################################################################################
#intoPopulation: putting data into existing population object (using intoPopulationSub.internal)
# 
# population - object of class population, data should be put into
# dataObject - matrix of data to be put into ril object
# dataType - what kind of data dataObject contains:
# 	-  founders - founders phenotypic
# 	-  offspring$phenotypes - offspring phenotypic
# 	-  offspring$genotypes - offspring genotype
# 	-  maps$genetic - genetic map 
# 	-  maps$physical - physical map
#
############################################################################################################
intoPopulation <- function(population, dataObject, dataType=c("founders","offspring$phenotypes","offspring$genotypes","maps$genetic","maps$physical"),verbose=FALSE,debugMode=0){
	### checks
	if(is.null(population)) stop("No population object provided!\n")
	inListCheck.internal(dataType,"dataType",c("founders","offspring$phenotypes","offspring$genotypes","maps$genetic","maps$physical"))
	if(verbose && debugMode==1) cat("intoPopulation starting without errors in checkpoints.\n")
	s <- proc.time()
	if(length(dataType)==1) {
		if(class(dataObject)!="list") stop("Multiple dataObjects should be provided as list.\n")
		if(length(dataObject)!=length(dataType)) stop("Support dataType for every element of dataObject.\n")
		if(length(dataType)!=length(unique(dataType))) stop("Every element of dataType must be unique!\n")
		for(i in 1:length(dataObject)){
			population <- intoPopulationSub.internal(population,dataObject[[i]],dataType[i], verbose, debugMode)
		}
	}
	else if(length(dataType)>1){
		population <- intoPopulationSub.internal(population,dataObject,dataType, verbose, debugMode)
	}

	if(is.null(population)) stop("No data provided!\n")
	class(population) <- "population"
	e <- proc.time()
	cat("intoPopulation for",dataType,"done in:",(e-s)[3],"seconds.\n")
	invisible(population)
}

############################################################################################################
#intoPopulationSub.internal: subfunction of intoPopulation, adding data to population object
# 
# population - object of class population, data should be put into
# dataObject - matrix of data to be put into ril object
# dataType - what kind of data dataObject contains:
# 	-  founders - founders phenotype
# 	-  offspring$phenotypes - offspring phenotype
# 	-  offspring$genotypes - offspring genotype
# 	-  maps$genetic - genetic map 
# 	-  maps$physical - physical map
#
############################################################################################################
intoPopulationSub.internal <- function(population, dataObject, dataType=c("founders","offspring$phenotypes","offspring$genotypes","maps$genetic","maps$physical"),verbose=FALSE,debugMode=0){
	if(verbose && debugMode==1) cat("intoPopulationSub.internal starting.\n")
	s <- proc.time()
	if(!(is.null(dataObject))&&!is.null(dim(dataObject))){
		### initialization
		columns <- NULL
		rows <- NULL
		
		### checking whether columns are numeric/convertable to numeric
		for(column in 1:ncol(dataObject)){
			if(!(numericCheck.internal(dataObject[,column],allow.na=TRUE))){
				columns <- c(columns,column)
			}
		}
		
		### removing faulty columns
		if(!(is.null(columns))){
			cat("Following  columns are not numeric and cannot be converted into numeric:",columns," so will be removed.\n")
			dataObject <- dataObject[,-columns]
		}
		
		if(is.null(dim(dataObject))) stop("Not enough data to continue.\n")
		
		### checking whether rows are numeric/convertable to numeric
		for(row_ in 1:nrow(dataObject)){
			if(!(numericCheck.internal(dataObject[row_,],allow.na=TRUE))){
				rows <- c(rows,row_)
			}
		}
		
		### removing faulty rows
		if(!(is.null(rows))){
			cat("Following  rows are not numeric and cannot be converted into numeric:",rows," so will be removed.\n")
			dataObject <- dataObject[-rows,]
		}
		
		if(is.null(dim(dataObject))) stop("Not enough data to continue.\n")
		
		cur<- matrix(as.numeric(as.matrix(dataObject)),nrow(dataObject),ncol(dataObject))
		
		### keeping colnames
		if(!is.null(colnames(dataObject))){
			colnames(cur) <- colnames(dataObject)
		}else{
			colnames(cur) <- 1:ncol(cur)
		}
		
		### keeping rownames
		if(!is.null(rownames(dataObject))){
			rownames(cur) <- rownames(dataObject)
		}else{
			rownames(cur) <- 1:nrow(cur)
		}
		
		### adding data to population
		if(dataType=="founders"){
			population$founders$phenotypes <- cur
			population$parameters$intoRil$founders <- list("dataObject", dataType, selectedRows, selectedCols)
			names(population$parameters$intoRil$founders) <- c("dataObject", "dataType", "selectedRows", "selectedCols")
		}else if(dataType=="offspring$phenotypes"){
			population$offspring$phenotypes <- cur
			population$parameters$intoRil$offspring$phenotypes <- list("dataObject", dataType, selectedRows, selectedCols)
			names(population$parameters$intoRil$offspring$phenotypes) <- c("dataObject", "dataType", "selectedRows", "selectedCols")
		}else if(dataType=="offspring$genotypes"){
			population$offspring$genotypes$real <- cur
			population$parameters$intoRil$offspring$genotypes <- list("dataObject", dataType, selectedRows, selectedCols)
			names(population$parameters$intoRil$offspring$genotypes) <- c("dataObject", "dataType", "selectedRows", "selectedCols")
		}else if(dataType=="maps$genetic"){
			population$maps$genetic <- cur
			population$parameters$intoRil$maps$genetic <- list("dataObject", dataType, selectedRows, selectedCols)
			names(population$parameters$intoRil$maps$genetic) <- c("dataObject", "dataType", "selectedRows", "selectedCols")
		}else if(dataType=="maps$physical"){
			population$maps$physical <- cur
			population$parameters$intoRil$maps$physical <- list("dataObject", dataType, selectedRows, selectedCols)
			names(population$parameters$intoRil$maps$physical) <- c("dataObject", "dataType", "selectedRows", "selectedCols")
		}
		
	}else{
		stop("No data provided for ",dataType,"!\n")
	}
	e <- proc.time()
	cat("intoPopulation for",dataType,"done in:",(e-s)[3],"seconds.\n")
	invisible(population)
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
