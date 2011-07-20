############################################################################################################
#
# intoPopulation.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified July 2011
# last modified in version: 0.8.6
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
# Contains: createPopulation, intoPopulation, intoPopulationSub.internal, intoPopulationSubPheno.internal, intoPopulationSubGeno.internal, intoPopulationSubMap.internal   
#
############################################################################################################

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
createPopulation <- function(offspring_phenotypes, founders, founders_groups, offspring_genotypes, maps_genetic, maps_physical, no.warn=FALSE, verbose=FALSE,debugMode=0){
	if(verbose && debugMode==1) cat("createPopulation starting.\n")
	s <- proc.time()
	population <- NULL
	if(missing(offspring_phenotypes)){
		stop("No offspring phenotype data provided!\n")
	}else{
		population <- intoPopulationSub.internal(population, offspring_phenotypes, "offspring$phenotypes", verbose, debugMode)
	}
	if(missing(founders)){
		stop("No founders phenotype data provided!\n")
	}else{
		population <- intoPopulationSub.internal(population, founders, "founders", verbose, debugMode)
	}
	if(missing(founders_groups)){
		stop("No information about founders groups provided!\n")
	}else{
		population <- intoPopulationSub.internal(population, founders_groups, "founders$groups", verbose, debugMode)
	}
	if(missing(offspring_genotypes)){
		if(verbose && !(no.warn))cat("No offspring genotypic data provided. You can supply it later using intoPopulation.\n")
	}else{
		population <- intoPopulationSub.internal(population, offspring_genotypes, "offspring$genotypes", verbose, debugMode)
	}
	if(missing(maps_genetic)){
		if(verbose && !(no.warn))cat("No genotic map provided. You can supply it later using intoPopulation.\n")
	}else{
		population <- intoPopulationSub.internal(population, maps_genetic, "maps$genetic", verbose, debugMode)
	}
	if(missing(maps_physical)){
		if(verbose && !(no.warn))cat("No physical map provided.  You can supply it later using intoPopulation.\n")
	}else{
		population <- intoPopulationSub.internal(population, maps_physical, "maps$physical", verbose, debugMode)
	}
	if(is.null(population)) stop("No data provided!\n")
	class(population) <- "population"
	e <- proc.time()
	is.population(population)
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
intoPopulation <- function(population, dataObject, dataType=c("founders","offspring$phenotypes","founders$group","offspring$genotypes","maps$genetic","maps$physical"),verbose=FALSE,debugMode=0){
	### checks
	is.population(population)
	inListCheck.internal(dataType,"dataType",c("founders","offspring$phenotypes","offspring$genotypes","maps$genetic","maps$physical"))
	if(verbose && debugMode==1) cat("intoPopulation starting without errors in checkpoints.\n")
	if(length(dataType)>1) {
		if(class(dataObject)!="list") stop("Multiple dataObjects should be provided as list.\n")
		if(length(dataObject)!=length(dataType)) stop("Support dataType for every element of dataObject.\n")
		if(length(dataType)!=length(unique(dataType))) stop("Every element of dataType must be unique!\n")
		for(i in 1:length(dataObject)){
			population <- intoPopulationSub.internal(population,dataObject[[i]],dataType[i], verbose, debugMode)
		}
	}
	else if(length(dataType)==1){
		population <- intoPopulationSub.internal(population,dataObject,dataType, verbose, debugMode)
	}

	if(is.null(population)) stop("No data provided!\n")
	class(population) <- "population"
	invisible(population)
}

############################################################################################################
#intoPopulationSubPheno.internal: subfunction of intoPopulation, adding data to population object
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
intoPopulationSub.internal <- function(population, dataObject, dataType=c("founders","offspring$phenotypes","founders$groups","offspring$genotypes","maps$genetic","maps$physical"),verbose=FALSE,debugMode=0){
	if(dataType=="founders" || dataType=="offspring$phenotypes"){
		population <- intoPopulationSubPheno.internal(population,dataObject,dataType, verbose, debugMode)
	}else if(dataType=="offspring$genotypes"){
		population <- intoPopulationSubGeno.internal(population,dataObject, verbose, debugMode)
	}else if(dataType=="maps$genetic"||dataType=="maps$physical"){
		population <- intoPopulationSubMap.internal(population,dataObject, dataType, verbose, debugMode)
	}else if(dataType=="founders$groups"){
		population$founders$groups <- dataObject
	}
	invisible(population)
}

############################################################################################################
#intoPopulationSubPheno.internal: subfunction of intoPopulation, adding data to population object
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
intoPopulationSubPheno.internal <- function(population, dataObject, dataType=c("founders","offspring$phenotypes"),verbose=FALSE,debugMode=0){
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
			if(verbose) cat("Following  columns are not numeric and cannot be converted into numeric:",columns," so will be removed.\n")
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
			if(verbose)cat("Following  rows are not numeric and cannot be converted into numeric:",rows," so will be removed.\n")
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
			population$parameters$intoRil$founders <- list("population object", "data object", dataType)
			names(population$parameters$intoRil$founders) <- c("population", "dataObject", "dataType")
		}else if(dataType=="offspring$phenotypes"){
			population$offspring$phenotypes <- cur
			population$parameters$intoRil$offspring$phenotypes <- list("population object", "data object", dataType)
			names(population$parameters$intoRil$offspring$phenotypes) <- c("population", "dataObject", "dataType")
		}
	}else{
		stop("No data provided for ",dataType,"!\n")
	}
	e <- proc.time()
	if(verbose&&debugMode==2)cat("intoPopulation for",dataType,"done in:",(e-s)[3],"seconds.\n")
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
intoPopulationSubGeno.internal <- function(population, dataObject, verbose=FALSE,debugMode=0){
	if(verbose && debugMode==1) cat("intoPopulationSubGeno.internal starting.\n")
	s <- proc.time()
	if(!(is.null(dataObject))&&!is.null(dim(dataObject))){
		### initialization
		columns <- NULL
		rows <- NULL
		
		### checking whether columns are numeric/convertable to numeric
		for(column in 1:ncol(dataObject)){
			if(!(genotypeCheck.internal(dataObject[,column],allow.na=TRUE))){
				columns <- c(columns,column)
			}
		}
		
		### removing faulty columns
		if(!(is.null(columns))){
			if(verbose)cat("Following  columns are not numeric and cannot be converted into numeric:",columns," so will be removed.\n")
			dataObject <- dataObject[,-columns]
		}
		
		if(is.null(dim(dataObject))) stop("Not enough data to continue.\n")
		
		### checking whether rows are numeric/convertable to numeric
		for(row_ in 1:nrow(dataObject)){
			if(!(genotypeCheck.internal(dataObject[row_,],allow.na=TRUE))){
				rows <- c(rows,row_)
			}
		}
		
		### removing faulty rows
		if(!(is.null(rows))){
			if(verbose)cat("Following  rows are not numeric and cannot be converted into numeric:",rows," so will be removed.\n")
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
		population$offspring$genotypes$real <- cur
		
	}else{
		stop("No data provided for offspring$genotypes !\n")
	}
	e <- proc.time()
	if(verbose&&debugMode==2)cat("intoPopulation for offspring$genotypes done in:",(e-s)[3],"seconds.\n")
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
intoPopulationSubMap.internal <- function(population, dataObject, dataType=c("maps$genetic","maps$physical"),verbose=FALSE, debugMode=0){
	if(verbose && debugMode==1) cat("intoPopulationSub.internal starting.\n")
	s <- proc.time()
	if(!(is.null(dataObject))&&!(is.null(dim(dataObject)))&&class(dataObject)=="matrix"){
		if(ncol(dataObject)!=2) stop ("This is not a correct map object.\n")
		if(any(!(is.numeric(dataObject)))) stop ("This is not a correct map object.\n")
		### adding data to population
		if(dataType=="maps$genetic"){
			population$maps$genetic <- dataObject
			population$parameters$intoRil$maps$genetic <- list("population object", "data object", dataType)
			names(population$parameters$intoRil$maps$genetic) <- c("population", "dataObject", "dataType")
		}else if(dataType=="maps$physical"){
			population$maps$physical <- dataObject
			population$parameters$intoRil$maps$physical <- list("population object", "data object", dataType)
			names(population$parameters$intoRil$maps$physical) <- c("population", "dataObject", "dataType")
		}
		
	}else{
		stop("No data provided for ",dataType,"!\n")
	}
	e <- proc.time()
	if(verbose&&debugMode==2)cat("intoPopulation for",dataType,"done in:",(e-s)[3],"seconds.\n")
	invisible(population)
}