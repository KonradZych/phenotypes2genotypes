#
# add.to.population.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Nov, 2012
# first written Mar, 2011
# Contains: create.population, add.to.population, add.to.populationSub.internal
#           add.to.populationSubPheno.internal, add.to.populationSubGeno.internal
#           add.to.populationSubMap.internal   
#

#  create.population
#
# DESCRIPTION:
#  Creating an object of class population using data supplied by user
# PARAMETERS:
#   - offspring$phenotypes - matrix containing offspring phenotype data (have to be supported, if not - function
#   - quits with error
#   - founders - matrix containing founders phenotype data (optional)
#   - offspring$genotypes - matrix containing offspring genotype data (optional)
#   - maps$genetic - matrix containing genetic map (optional)
#   - maps$physical - matrix containing physical map (optional)
# OUTPUT:
#  An object of class population
#
create.population <- function(offspring_phenotypes, founders, founders_groups, offspring_genotypes, maps_genetic, maps_physical,
  populationType=c("riself", "f2", "bc", "risib"), no.warn=FALSE, verbose=FALSE, debugMode=0){
  if(verbose && debugMode==1) cat("create.population starting.\n")
  s <- proc.time()
  population <- NULL
  populationType <- match.arg(populationType)
  if(missing(offspring_phenotypes)){
    stop("No offspring phenotype data provided!\n")
  }else{
    population <- add.to.populationSub.internal(population, populationType, offspring_phenotypes, "offspring$phenotypes", verbose, debugMode)
  }
  if(missing(founders)){
    cat("No founders phenotype data provided, it will be simulated!\n")
    offsprings <- population$offspring$phenotypes
    half <-floor(ncol(offsprings)/2)
    founders <- t(apply(offsprings,1,function(x){c(mean(sort(x)[1:half]),mean(sort(x)[2:(half+1)]),
    mean(sort(x)[3:(half+2)]),mean(sort(x)[(half+1):ncol(offsprings)]),mean(sort(x)[(half):ncol(offsprings)-1]),
    mean(sort(x)[(half-1):ncol(offsprings)-2]))}))
    population$flags <- c(population$flags,"noParents")
    founders_groups <- c(0,0,0,1,1,1)
    population <- add.to.populationSub.internal(population, populationType, founders, "founders", verbose, debugMode)
  }else{
    n.childrenNotInParental <- sum(!(rownames(founders)%in%rownames(population$offspring$phenotypes)))
    if(n.childrenNotInParental==nrow(founders)){
      stop("No match between the row names in the founders and offspring.\n")
    }else if(n.childrenNotInParental!=0){
      warning(n.childrenNotInParental,"markers from founders file are not present in offspring data and will be removed.\n")
      founders <- founders[which((rownames(founders)%in%rownames(population$offspring$phenotypes))),]
    }
    population <- add.to.populationSub.internal(population, populationType,founders, "founders", verbose, debugMode)
  }
  if(missing(founders_groups)){
    stop("No information about founders groups provided!\n")
  }else{
    if(length(founders_groups)==ncol(population$founders$phenotypes)){
      population <- add.to.populationSub.internal(population, populationType,founders_groups, "founders$groups", verbose, debugMode)
    }else{
      stop("founders_group parameter should have length equall to number of columns in founders phenotype data!\n")
    }
  }
  if(missing(offspring_genotypes)){
    if(verbose && !(no.warn))cat("No offspring genotypic data provided. You can supply it later using add.to.population.\n")
  }else{
    population <- add.to.populationSub.internal(population, populationType, offspring_genotypes, "offspring$genotypes", verbose, debugMode)
  }
  if(missing(maps_genetic)){
    if(verbose && !(no.warn))cat("No genotic map provided. You can supply it later using add.to.population.\n")
  }else{
    population <- add.to.populationSub.internal(population, populationType, maps_genetic, "maps$genetic", verbose, debugMode)
  }
  if(missing(maps_physical)){
    if(verbose && !(no.warn))cat("No physical map provided.  You can supply it later using add.to.population.\n")
  }else{
    population <- add.to.populationSub.internal(population, populationType, maps_physical, "maps$physical", verbose, debugMode)
  }
  if(is.null(population)) stop("No data provided!\n")
  class(population) <- c("population", populationType)
  e <- proc.time()
  check.population(population)
  if(verbose){
    if(debugMode==2){
      cat("create.population done in:",(e-s)[3],"seconds.\n")
    }else{
      cat("create.population finished\n")
    }
  }
  invisible(population)
}

############################################################################################################
#                  *** add.to.population ***
#
# DESCRIPTION:
#  putting data into existing population object (using add.to.populationSub.internal)
# 
# PARAMETERS:
#   population - object of class population, data should be put into
#   dataObject - matrix of data to be put into the population object
#   dataType - what kind of data dataObject contains:
#     -  founders - founders phenotypic
#     -  offspring$phenotypes - offspring phenotypic
#     -  offspring$genotypes - offspring genotype
#     -  maps$genetic - genetic map 
#     -  maps$physical - physical map
#   verbose - be verbose
#   debugMode - 1: print our checks, 2: print additional time information
#
# OUTPUT:
#  an object of class population
#
############################################################################################################
add.to.population <- function(population, dataObject, dataType=c("founders","offspring$phenotypes","founders$group","offspring$genotypes","maps$genetic","maps$physical"),verbose=FALSE,debugMode=0){
  ### checks
  check.population(population)
  populationType <- class(population)[2]
  dataType <- match.arg(dataType)
  if(verbose && debugMode==1) cat("add.to.population starting without errors in checkpoints.\n")
  if(length(dataType)>1){
    if(class(dataObject)!="list") stop("Multiple dataObjects should be provided as list.\n")
    if(length(dataObject)!=length(dataType)) stop("Support dataType for every element of dataObject.\n")
    if(length(dataType)!=length(unique(dataType))) stop("Every element of dataType must be unique!\n")
    for(i in 1:length(dataObject)){
      population <- add.to.populationSub.internal(population,populationType,dataObject[[i]],dataType[i], verbose, debugMode)
    }
  }else if(length(dataType)==1){
    population <- add.to.populationSub.internal(population,populationType,dataObject,dataType, verbose, debugMode)
  }

  if(is.null(population)) stop("No data provided!\n")
  invisible(population)
}

############################################################################################################
#                  *** add.to.populationSub.internal ***
#
# DESCRIPTION:
#  subfunction of add.to.population, using subfunctions to add a single data object to the object of class
#  population
# 
# PARAMETERS:
#   population - object of class population, data should be put into
#   dataObject - matrix of data to be put into the population object
#   dataType - what kind of data dataObject contains:
#     -  founders - founders phenotype
#     -  offspring$phenotypes - offspring phenotype
#     -  offspring$genotypes - offspring genotype
#     -  maps$genetic - genetic map 
#     -  maps$physical - physical map
#   verbose - be verbose
#   debugMode - 1: print our checks, 2: print additional time information
#
# OUTPUT:
#  an object of class population
#
############################################################################################################
add.to.populationSub.internal <- function(population, populationType=c("riself", "f2", "bc", "risib"), dataObject, dataType=c("founders","offspring$phenotypes","founders$groups","offspring$genotypes","maps$genetic","maps$physical"),verbose=FALSE,debugMode=0){
  dataType <- match.arg(dataType)
  populationType <- match.arg(populationType)
  if(missing(dataObject)) stop("dataObject is missing\n")
  if(dataType!="founders$groups"){
    if(class(dataObject)=="data.frame"){
      dataObject <- as.matrix(dataObject)
    }else if(class(dataObject)!="matrix"){
      stop("dataObject should be either a matrix or a date frame")
    }
  }
  if(dataType=="founders" || dataType=="offspring$phenotypes"){
    population <- add.to.populationSubPheno.internal(population,dataObject,dataType, verbose, debugMode)
  }else if(dataType=="offspring$genotypes"){
    if(!(is.null(dataObject))&&!is.null(dim(dataObject))){  
      population$offspring$genotypes$real <- add.to.populationSubGeno.internal(population,dataObject,populationType,verbose)
    }else{
      stop("No data provided for offspring$genotypes !\n")
    }
  }else if(dataType=="maps$genetic"||dataType=="maps$physical"){
    population <- add.to.populationSubMap.internal(population,dataObject, dataType, verbose, debugMode)
  }else if(dataType=="founders$groups"){
    population$founders$groups <- dataObject
  }
  invisible(population)
}

############################################################################################################
#                  *** add.to.populationSubPheno.internal ***
#
# DESCRIPTION:
#  subfunction of add.to.populationSub.internal, adding a single phenotype object to the object of class 
#  population
# 
# PARAMETERS:
#   population - object of class population, data should be put into
#   dataObject - matrix of data to be put into the population object
#   dataType - what kind of data dataObject contains:
#     -  founders - founders phenotype
#     -  offspring$phenotypes - offspring phenotype
#   verbose - be verbose
#   debugMode - 1: print our checks, 2: print additional time information
#
# OUTPUT:
#  an object of class population
#
############################################################################################################
add.to.populationSubPheno.internal <- function(population, dataObject, dataType=c("founders","offspring$phenotypes"),verbose=FALSE,debugMode=0){
  if(verbose && debugMode==1) cat("add.to.populationSub.internal starting.\n")
  s <- proc.time()
  if(!(is.null(dataObject))&&!is.null(dim(dataObject))){
    #Check whether rows are numeric/convertable to numeric
    rows <- unlist(lapply(c(1:nrow(dataObject)),add.to.populationSubPhenoSub.internal,dataObject,verbose))    
    #Removing faulty rows
    if(!(is.null(rows))){
      if(verbose)cat("Following  rows are not numeric and cannot be converted into numeric:",rows," so will be removed.\n")
      dataObject <- dataObject[-rows,]
    }
    
    if(is.null(dim(dataObject))) stop("Not enough data to continue.\n")
    
    cur<- apply(as.matrix(dataObject),c(1,2),as.numeric)
    
    #Keeping colnames
    if(!is.null(colnames(dataObject))){
      colnames(cur) <- colnames(dataObject)
    }else{
      colnames(cur) <- 1:ncol(cur)
    }
    
    #Keeping rownames
    if(!is.null(rownames(dataObject))){
      rownames(cur) <- rownames(dataObject)
    }else{
      rownames(cur) <- 1:nrow(cur)
    }
    
    #Adding data to population
    if(dataType=="founders"){
      population$founders$phenotypes <- cur
    }else if(dataType=="offspring$phenotypes"){
      population$offspring$phenotypes <- cur
    }
  }else{
    stop("No data provided for ",dataType,"!\n")
  }
  e <- proc.time()
  if(verbose&&debugMode==2)cat("add.to.population for",dataType,"done in:",(e-s)[3],"seconds.\n")
  invisible(population)
}

############################################################################################################
#                  *** add.to.populationSubPhenoSub.internal ***
#
# DESCRIPTION:
#  subfunction of add.to.populationSubPheno.internal, checking if a single row of dataObject is formatted
#  correctly
# 
# PARAMETERS:  
#  curRow - number of row currently being checked  
#   dataObject - matrix of data to be put into the population object
#   verbose - be verbose (show information about the progress)
#
# OUTPUT:
#  number or NULL
#
############################################################################################################
add.to.populationSubPhenoSub.internal <- function(curRow,dataObject,verbose){
  if(verbose&&curRow%%1000==0) cat("Processing row:",curRow,"\n")
  if(!(numericCheck.internal(dataObject[curRow,],allow.na=TRUE))){
    return(curRow)
  }else{
    return(NULL)
  }
}

############################################################################################################
#                  *** add.to.populationSubGenoSub.internal ***
#
# DESCRIPTION:
#  subfunction of add.to.populationSubGeno.internal, checking if a single row of dataObject is formatted
#  correctly
# 
# PARAMETERS:  
#  curRow - number of row currently being checked  
#   dataObject - matrix of data to be put into the population object
#   verbose - be verbose (show information about the progress)
#
# OUTPUT:
#  number or NULL
#
############################################################################################################
add.to.populationSubGeno.internal <- function(population,dataObject,populationType=c("riself", "f2", "bc", "risib"),verbose=FALSE){
  populationType <- match.arg(populationType)
    #Checking whether rows are numeric/convertable to numeric
      if(populationType == "f2"){
        genotypes <- c(1:5)
      }else{
        genotypes <- c(1,2)
      }
      rows <- unlist(lapply(c(1:nrow(dataObject)),add.to.populationSubGenoSub.internal,dataObject,genotypes,verbose))
      #Removes faulty rows
      if(!(is.null(rows))){
        if(verbose)cat("Following  rows are not numeric and cannot be converted into numeric:",rows," so will be removed.\n")
        dataObject <- dataObject[-rows,]
      }
    
      if(is.null(dim(dataObject))) stop("Not enough data to continue.\n")
    
      cur<- matrix(as.numeric(as.matrix(dataObject)),nrow(dataObject),ncol(dataObject))
    
      #Keep colnames
      if(!is.null(colnames(dataObject))){
        colnames(cur) <- colnames(dataObject)
      }else{
        colnames(cur) <- 1:ncol(cur)
      }
    
      #Keep rownames
      if(!is.null(rownames(dataObject))){
        rownames(cur) <- rownames(dataObject)
      }else{
        rownames(cur) <- 1:nrow(cur)
      }
      if(verbose){
        for(x in as.numeric(names(table(cur)))){
          cat(x,": ",round(sum(cur==x,na.rm=T)/length(cur)*100,2),"%\n",sep="")
        }
        cat("NA: ",round(sum(is.na(cur))/length(cur)*100,2),"%\n",sep="")
      }
      #Adding data to population
      return(cur)
}

add.to.populationSubGenoSub.internal <- function(curRow,dataObject,genotypes,verbose=FALSE){
  if(verbose&&curRow%%1000==0) cat("Processing row:",curRow,"\n")
  if(!(genotypeCheck.internal(dataObject[curRow,],genotypes,allow.na=TRUE))){
    return(curRow)
  }else{
    return(NULL)
  }
}

############################################################################################################
#                  *** add.to.populationSubMap.internal ***
#
# DESCRIPTION:
#  subfunction of add.to.populationSub.internal, adding a single map object to the object of class population
# 
# PARAMETERS:  
#  population - object of class population, data should be put into
#   dataObject - matrix of data to be put into ril object
#   dataType - what kind of data dataObject contains:
#     -  maps$genetic - genetic map 
#     -  maps$physical - physical map
#   verbose - be verbose
#   debugMode - 1: print our checks, 2: print additional time information
#
# OUTPUT:
#  an object of class population
#
############################################################################################################
add.to.populationSubMap.internal <- function(population, dataObject, dataType=c("maps$genetic","maps$physical"),verbose=FALSE, debugMode=0){
  if(verbose && debugMode==1) cat("add.to.populationSub.internal starting.\n")
  s <- proc.time()
  if(!(is.null(dataObject))&&!(is.null(dim(dataObject)))&&class(dataObject)=="matrix"){
    if(dataType=="maps$genetic" && ncol(dataObject)!=2) cat ("This is not a correct map object.\n")
    if(dataType=="maps$physical" && ncol(dataObject)!=2 && ncol(dataObject)!=3) cat ("This is not a correct map object.\n")
    if(any(!(is.numeric(dataObject)))) stop ("This is not a correct map object.\n")
    ### adding data to population
    if(dataType=="maps$genetic"){
      population$maps$genetic <- dataObject
    }else if(dataType=="maps$physical"){
      if(ncol(dataObject) == 3){
        population$maps$physical <- dataObject
        colnames(population$maps$physical) <- c("Chr","Start","End")
      }else{
        population$maps$physical <- cbind(dataObject,dataObject[,2])
        colnames(population$maps$physical) <- c("Chr","Start","End")
      }
    }
  }else{
    stop("No data provided for ",dataType,"!\n")
  }
  e <- proc.time()
  if(verbose&&debugMode==2)cat("add.to.population for",dataType,"done in:",(e-s)[3],"seconds.\n")
  invisible(population)
}
