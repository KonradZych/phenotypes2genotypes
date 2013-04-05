#
# readFiles.r
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Nov, 2012
# first written Nov, 2011
# Contains: read.population mapMarkers.internal, correctRowLoc.internal
#           probesLocation.internal, orrectRowGff.internal
#

#  read.population
#
# DESCRIPTION:
#  Reads geno/phenotypic files into R environment into special object.
# PARAMETERS:
#   - offspring - Core used to specify names of children phenotypic ("offspring_phenotypes.txt") and genotypic ("offspring_genotypes.txt") files.
#   - founders - Core used to specify names of founders phenotypic ("founders_phenotypes.txt") file.
#   - map - Core used to specify names of genetic ("map_genetic.txt") and physical ("map_physical.txt") map files.
#   - founders_groups - specify founders groups
#   - verbose - Be verbose
#   - debugMode - 1: Print our checks, 2: print additional time information
# OUTPUT:
#   An object of class population 
#
read.population <- function(offspring = "offspring", founders = "founders", map = "map", foundersGroups, populationType = c("riself", "f2", "bc", "risib"),
  readMode = c("normal","HT"), verbose = FALSE, debugMode = 0, ...){

  ### checks
  populationType <- match.arg(populationType)
  readMode <- match.arg(readMode)
  
  
  ### file names
  fileFoundersPheno  <- paste(founders,"_phenotypes.txt",sep="")
  fileOffspringPheno <- paste(offspring,"_phenotypes.txt",sep="")
  fileOffspringGeno  <- paste(offspring,"_genotypes.txt",sep="")
  fileAnnotations    <- paste(offspring,"_annotations.txt",sep="")
  fileMapPhys        <- paste(map,"_physical.txt",sep="")
  fileMapGen         <- paste(map,"_genetic.txt",sep="")

  ### initializing  
  s <- proc.time()
  if(verbose && debugMode==1) cat("read.population starting.\n")
  population <- NULL
  
  ### offspring phenotypic file
  if(!file.exists(fileOffspringPheno)){
    stop("No phenotype file for offspring: ",fileOffspringPheno," this file is essential, you have to provide it\n")
  }else{
    if(verbose)  cat("File:",fileOffspringPheno,"found and will be processed.\n")
    if(readMode == "normal"){
      ### TODO: this should be using readSingleFile
      offspringPhenotypes <- read.table(fileOffspringPheno,sep="\t", row.names=1, header=TRUE)
      population <- add.to.populationSub.internal(population,offspringPhenotypes,"offspring$phenotypes",populationType=populationType)
    }else{
      population$offspring$phenotypes <- fileOffspringPheno
    }
  }

  ### founders phenotypic file
  if(!file.exists(fileFoundersPheno)){
    ### simulate data if there is no file
    if(verbose)cat("No phenotype file for founders: ",fileFoundersPheno,". Founder phenotypes will be simulated.\n")
    if(readMode == "normal"){
      population <- simulateParentalPhentypes(population, population$offspring$phenotypes, populationType)
    }
    ### if the mode is HT, we don't simulate founders but just judge the variance in the offspring while converting phenotypes to genotypes
  }else{
    ### read the file if present
    if(verbose)  cat("File:",fileFoundersPheno,"found and will be processed.\n")
    ### check if there is an information about the founders groups
    if(missing(foundersGroups)) stop("No information about founders groups provided.\n")
    ### founders groups should be a sequence of 0s and 1s
    if(any(foundersGroups!=0 && foundersGroups!=1)) stop("Founders groups attribute is incorrect.\n")
    if(readMode == "normal"){
      population <- readSingleFile(population, fileFoundersPheno, "founders", verbose=verbose, header=TRUE)
      if(length(foundersGroups)!=ncol(population$offspring$phenotypes)) stop("Founders groups attribute is incorrect.\n")
    }else{
      population <- readFoundersAndTtest(fileFoundersPheno, founders_groups, populationType, threshold, sliceSize, transformations, verbose)
    }
  }
  
  class(population) <- c("population",populationType)
  
  ### annotations file
  population <- readSingleFile(population, fileAnnotations, "annotations", verbose=verbose, header=TRUE)

  ### offspring genotypic file
  population <- readSingleFile(population, fileOffspringGeno, "offspring$genotypes", verbose=verbose, header=TRUE)
  
  ### physical map
  population <- readSingleFile(population, fileMapPhys, "maps$physical", verbose=verbose, header=FALSE)

  ### genetic map
  population <- readSingleFile(population, fileMapGen, "maps$genetic", verbose=verbose, header=FALSE)

  ### TO BE REMOVED when generate.biomarkers is corrected
  population$sliceSize <- 5000
  #**********FINALIZING FUNCTION*************
  e <- proc.time()
  if(verbose && debugMode==2) cat("read.population finished after",(e-s)[3],"seconds.\n")
  invisible(population)
}

#  readSingleFile
#
# DESCRIPTION:
#  Reads single geno/phenotypic file into R environment into special object.
# PARAMETERS:
#   - population - an object of class population
#   - filename - name of the file that will be processed
#   - fileType - type of the file that will be processed (information about what type of data it contains)
# OUTPUT:
#   An object of class population 
#
readSingleFile <- function(population, filename, fileType, verbose=FALSE, ...){
  if(file.exists(filename)){
    if(verbose) cat("File:",filename,"found and will be processed.\n")
    dataRead <- read.table(filename,sep="\t", row.names=1, ...)
    population <- add.to.population(population, dataRead, fileType)
    doCleanUp.internal()
  }else{
    if(verbose) cat("File:",filename,"not found.\n")
  }
  invisible(population)
}

TwoClassTtest <- function(foundersline,foundersGroups,threshold){
  p1  <- as.numeric(foundersline[which(foundersGroups == 0) + 1])  #TODO: Why add 1, Why dont we get a out-of-bounds error ?? -> first element of the line is a name of the row
  p2  <- as.numeric(foundersline[which(foundersGroups == 1) + 1])  #TODO: Why add 1, Why dont we get a out-of-bounds error ??
  res <- t.test(p1, p2)
  if(res$p.value < threshold) return(foundersline)
  return(NULL)
}


simulateParentalPhentypes <- function(population, offspringPhenotypes, populationType){
  cat("No founders phenotype data provided, it will be simulated!\n")
  half     <- floor(ncol(offspringPhenotypes)/2)
  end      <- ncol(offspringPhenotypes)
  founders <- t(apply(offspringPhenotypes, 1, function(x){
    x <- sort(x)
    c( mean(x[1:half],na.rm=TRUE), mean(x[2:(half+1)],na.rm=TRUE), mean(x[3:(half+2)],na.rm=TRUE),
       mean(x[(half):end],na.rm=TRUE), mean(x[(half-1):(end-1)],na.rm=TRUE), mean(x[(half-2):end-2],na.rm=TRUE))
  }))
  population$flags <- c(population$flags,"noParents")
  population <- add.to.populationSub.internal(population, founders, "founders", populationType=populationType)
  population$founders$groups <- c(0,0,0,1,1,1)
  return(population)
}

byLineTtest <- function(filename, foundersGroups, threshold, sliceSize, transformations, verbose=FALSE){
  s <- proc.time()
  analysedFile <- file(filename,"r")
  if(verbose) cat("Analysing file:",filename,"\n")
  analysedLines <- NULL
  result        <- NULL
  count         <- 0
  n.columns     <- length(foundersGroups) + 1 # Why +1
  header        <- unlist(strsplit(readLines(analysedFile, n = 1), "\t"))
  if(length(header) != length(foundersGroups)){
    if(length(header) != (length(foundersGroups) + 1)){
      stop("Incorrect founders_groups parameter. Should have length: ",length(header)," instead of: ",length(founders_groups),"\n")
    }else{
      header <- header[-1]
    }
  }
  analysedLines <- readLines(analysedFile, n=sliceSize) # Assume we have at least 1 line
  while(!((length(analysedLines) == 0) && (typeof(analysedLines) == "character"))){
    readData  <- t(matrix(unlist(strsplit(analysedLines,"\t")), n.columns, length(analysedLines)))
    nReadData <- apply(readData[,-1], c(1,2), as.numeric)
    dataUsed  <- cbind(nReadData[,1], t(transformation(nReadData, transformations)[[length(transformations)]]))

    selected  <- rbind(unlist(apply(dataUsed,1,TwoClassTtest, foundersGroups, threshold)))
    selected  <- t(matrix(selected, n.columns, length(selected)/n.columns))
    oldnames  <- rownames(result)
    result    <- rbind(result, apply(selected[,-1],c(1,2),as.numeric))  # Isn't this Nread Data ????
    rownames(result) <- c(oldnames, selected[,1])
    count <- count + sliceSize
    if(verbose) cat("Processed: ",count,"lines\n")
    analysedLines <- readLines(analysedFile, n=sliceSize)
  }
  close(analysedFile)
  if((is.null(result))) stop("No parental data read in!\n")
  colnames(result) <- header
  invisible(result)
}

loadByLine <- function(filename,pattern,sliceSize,transformations,verbose=FALSE){
  s <- proc.time()
  analysedFile  <- file(filename,"r")
  result        <- NULL
  count         <- 0
  header        <- unlist(strsplit(readLines(analysedFile,n=1),"\t"))
  n.columns     <- length(header)+1
  analysedLines <- readLines(con=analysedFile, n=sliceSize)
  while(!((length(analysedLines) == 0) && (typeof(analysedLines) == "character"))){
    readData <- t(matrix(unlist(strsplit(analysedLines,"\t")), n.columns, length(analysedLines)))
    selected <- rbind(unlist(apply(readData, 1, function(x){
      if(x[1] %in% pattern) return(x)
      return(NULL)
    }, pattern)))
    selected <- t(matrix(selected,n.columns,length(selected)/n.columns))
    result   <- rbind(result, selected)
    count    <- count + sliceSize
    if(verbose) cat("Processed:", count, "lines\n")
    analysedLines <- readLines(con=analysedFile, n=sliceSize)
  }
  close(analysedFile)
  rownames(result) <- result[,1]
  result <- apply(result[,-1], c(1,2), as.numeric)
  if(transformations!="nothing") result <- transformation(result,transformations)[[length(transformations)]]
  colnames(result) <- header
  invisible(result)
}

readFoundersAndTtest <- function(fileFoundersPheno, founderGroups, populationType, threshold, sliceSize, transformations, verbose){
  if(file.exists(fileFoundersPheno)){
    population <- NULL
    if(verbose) cat("Found phenotypic file for founders:", fileFoundersPheno,"and will store  it in population$founders$phenotypes\n")
    foundersPhenotypes <- byLineTtest(fileFoundersPheno, founderGroups, threshold, sliceSize, transformations, verbose)
    population <- add.to.populationSub.internal(population, foundersPhenotypes, "founders", populationType=populationType)
    population$founders$groups <- founderGroups
    doCleanUp.internal()
    invisible(return(population))
  }else{
    cat("No phenotype file for founders at location:", fileFoundersPheno,"\n")
    #TODO: SIMULATE DATA WHEN NOT AVAILABLE !!!!!
  }
}

############################################################################################################
#                  *** mapMarkers.internal ***
#
# DESCRIPTION:
#  removes from matrix1 cols or rows, which are not present in second (coparing using col/rownames)
# 
# PARAMETERS:
#   expressionMatrix1, expressionMatrix2 - matrices with data of any type
#   mapMode - 1 - map rows, 2 - map cols
#   verbose - Be verbose
#   debugMode - 1: Print our checks, 2: print additional time information
#
# OUTPUT:
#  object of class population 
#
############################################################################################################
mapMarkers.internal <- function(expressionMatrix1, expressionMatrix2, mapMode=2, verbose=FALSE, debugMode=0){
  if(mapMode==1) {
    nrRows <- nrow(expressionMatrix1)
    ### warnings when names are mismatching
    if(verbose && debugMode==2)if(nrRows!=nrow(expressionMatrix2)){
      cat("Following markers will be removed:\n")
      cat(paste(rownames(expressionMatrix1)[which(!(rownames(expressionMatrix1) %in% rownames(expressionMatrix2)))],"\n"))
    }
    ### mapping itself
    expressionMatrix1 <- expressionMatrix1[which(rownames(expressionMatrix1) %in% rownames(expressionMatrix2)),]
    if(verbose) cat("Because of names mismatch,",nrRows-nrow(expressionMatrix1),"markers were removed, run function with verbose=T debugMode=2 to print their names out.\n")
  }
  else if(mapMode==2){
    nrCols <- ncol(expressionMatrix1)
    ### warnings when names are mismatching
    if(verbose && debugMode==2)if(nrCols!=ncol(expressionMatrix2)){
      cat("Following individuals will be removed:\n")
      paste(colnames(expressionMatrix1)[which(!(colnames(expressionMatrix1) %in% colnames(expressionMatrix2)))],"\n")
    }
    ### mapping itself
    expressionMatrix1 <- expressionMatrix1[,which(colnames(expressionMatrix1) %in% colnames(expressionMatrix2))]
    if(verbose) cat("Because of names mismatch,",nrCols-ncol(expressionMatrix1),"individuals were removed, run function with verbose=T debugMode=2 to print their names out.\n")
  }
  invisible(expressionMatrix1)
}
