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
read.population <- function(offspring="offspring",founders="founders",map="maps",founders_groups,populationType=c("riself", "f2", "bc", "risib"),
  read_mode=c("normal","HT"), verbose=FALSE,debugMode=0,...){
  s <- proc.time()
  if(verbose && debugMode==1) cat("read.population starting.\n")
  population <- NULL
  populationType <- match.arg(populationType)
  read_mode <- match.arg(read_mode)
  if(read_mode=="normal"){
    population <- read.populationNormal.internal(offspring, founders, map, founders_groups, populationType, verbose, debugMode)
  }else{
    population <- read.populationHT.internal(offspring, founders, map, founders_groups, populationType, verbose, debugMode,...)
  }
  #**********FINALIZING FUNCTION*************
  e <- proc.time()
  if(verbose) cat("read.population done in",(e-s)[3],"seconds.\n")
  class(population) <- c("population", populationType)
  doCleanUp.internal()
  invisible(population)
}

read.populationNormal.internal <- function(offspring,founders,map,founders_groups,populationType=c("riself", "f2", "bc", "risib"), verbose=FALSE,debugMode=0){
  #**********READING CHILDREN PHENOTYPIC DATA*************
  filename <- paste(offspring,"_phenotypes.txt",sep="")
  if(file.exists(filename)){
    if(verbose) cat("Found phenotypic file for offspring:",filename,"and will store  it in population$offspring$phenotypes\n")
    offspring_phenotypes <- read.table(filename,sep="\t",header=TRUE)
    offspring_phenotypes <- as.matrix(offspring_phenotypes)
    population <- NULL
    population <- add.to.populationSub.internal(population,offspring_phenotypes,"offspring$phenotypes",populationType=populationType)
    doCleanUp.internal()
  }else{
    stop("No phenotype file for offspring: ",filename," this file is essentiall, you have to provide it\n")
  }
  
  #**********READING PARENTAL PHENOTYPIC DATA*************
  filename <- paste(founders,"_phenotypes.txt",sep="")
  if(missing(founders)){
    cat("No founders phenotype data provided, it will be simulated!\n")
    offsprings <- population$offspring$phenotypes
    half <-floor(ncol(offsprings)/2)
    founders <- t(apply(offsprings,1,function(x){c(mean(sort(x)[1:half]),mean(sort(x)[2:(half+1)]),
    mean(sort(x)[3:(half+2)]),mean(sort(x)[(half+1):ncol(offsprings)]),mean(sort(x)[(half):ncol(offsprings)-1]),
    mean(sort(x)[(half-1):ncol(offsprings)-2]))}))
    population$flags <- c(population$flags,"noParents")
    population <- add.to.populationSub.internal(population, founders, "founders",populationType=populationType)
    founders_groups <- c(0,0,0,1,1,1)
  }else if(file.exists(filename)){
    if(verbose) cat("Found phenotypic file for founders:",filename,"and will store it in population$founders$phenotypes\n")
    founders <- read.table(filename,sep="\t",header=TRUE)
    founders <- as.matrix(founders)
    population <- add.to.populationSub.internal(population, founders, "founders",populationType=populationType)
    #removing from founders probes that are not in children:
    population$founders$phenotypes <- mapMarkers.internal(population$founders$phenotypes,population$offspring$phenotypes, mapMode=1, verbose=verbose)
    doCleanUp.internal()
  }else{
    stop("No phenotype file for founders at location: ",filename," \n")
  }
  
  if(missing(founders_groups)){ 
    stop("Specify founders_groups!\n")
  }
  #**********FOUNDERS GROUPS*************
  if(length(founders_groups)!=ncol(population$founders$phenotypes)){
    stop("Length of founders_groups should be equal to the number of columns in founder phenotype file.")
  }else if((sum(founders_groups==1)+sum(founders_groups==0))!=length(founders_groups)){
    stop("founders_groups should contain only 0s and 1s.")
  }else{
    population$founders$groups <- founders_groups
  }
  
  #**********READING CHILDREN GENOTYPIC DATA*************
  filename <- paste(offspring,"_genotypes.txt",sep="")
  if(file.exists(filename)){
    if(verbose) cat("Found genotypic file for offspring:",filename,"and will store  it in population$offspring$genotypes$real\n")
    offspring_genotypes <- read.table(filename,sep="\t",header=TRUE)
    offspring_genotypes <- as.matrix(offspring_genotypes)
    population <- add.to.population(population, offspring_genotypes, "offspring$genotypes",populationType=populationType)
    doCleanUp.internal()
  }else{
    if(verbose)cat("No genotypic file for offspring:",filename,"genotypic data for offspring will be simulated\n")
  }
  
  #**********READING GENETIC MAP*************
  filename <- paste(map,"_genetic.txt",sep="")
  if(file.exists(filename)){
    if(verbose) cat("Found genotypic file for offspring:",filename,"and will store  it in population$maps$genetic\n")
    maps_genetic <- read.table(filename,sep="\t",row.names=1,header=FALSE)
    maps_genetic <- as.matrix(maps_genetic)
    population <- add.to.population(population, maps_genetic, "maps$genetic",populationType=populationType)
    doCleanUp.internal()
  }else{
    if(verbose)cat("No genetic map file:",filename,".\n")
  }
  
  #**********READING PHYSICAL MAP*************
  filename <- paste(map,"_physical.txt",sep="")
  if(file.exists(filename)){
    if(verbose) cat("Found genotypic file for offspring:",filename,"and will store  it in population$maps$physical\n")
    physical <-  read.table(filename,sep="\t",row.names=1,header=FALSE)
    physical <- as.matrix(physical)
    population <- add.to.population(population, physical, "maps$physical",populationType=populationType)
    doCleanUp.internal()
  }else{
    if(verbose)cat("No physical map file:",filename,".\n")
  }
  invisible(population)
}

t.test.line.internal<- function(foundersline,founders_groups,threshold){
  parentA <- as.numeric(foundersline[which(founders_groups==0)+1])
  parentB <- as.numeric(foundersline[which(founders_groups==1)+1])
  res <- t.test(parentA,parentB)
  if(res$p.value < threshold){
    return(foundersline)
  }else{
    return(NULL)
  }
}

t.test.linebyline.internal <- function(filename,founders_groups,threshold,sliceSize,verbose=FALSE){
  s <- proc.time()
  analysedFile <- file(filename,"r")
  if(verbose)cat("Analysing file:",filename,"\n")
  analysedLines <- NULL
  result <- NULL
  count <- 0
  n.columns <- length(founders_groups)+1
  header <- unlist(strsplit(readLines(analysedFile,n=1),"\t"))
  if(length(header)!=length(founders_groups)){
    if(length(header)!=length(founders_groups)+1){
      stop("Incorrect founders_groups parameter. Should have length: ",length(header)," instead of: ",length(founders_groups),"\n")
    }else{
      header <- header[-1]
    }
  }
  while(!((length(analysedLines) == 0) && (typeof(analysedLines) == "character"))){
    analysedLines<-readLines(analysedFile,n=sliceSize)
    if(!((length(analysedLines) == 0) && (typeof(analysedLines) == "character"))){
      aa <- t(matrix(unlist(strsplit(analysedLines,"\t")),n.columns,length(analysedLines)))
      analysedLines <- NULL
      gc()
      selected <- rbind(unlist(apply(aa,1,t.test.line.internal,founders_groups,threshold)))
      selected <- t(matrix(selected,n.columns,length(selected)/n.columns))
      result <- rbind(result,selected)
      count <- count + sliceSize
      if(verbose) cat("Processed: ",count,"lines\n")
    }
  }
  if((is.null(result))){
    stop("No parental data read in!\n")
  }
  rownames(result) <- result[,1]
  result <- apply(result[,-1],c(1,2),as.numeric)
  colnames(result) <- header
  invisible(result)
}

load.line.internal <- function(curLine, pattern){
  if(curLine[1] %in% pattern){
    return(curLine)
  }
  return(NULL)
}

load.linebyline.internal <- function(filename,pattern,sliceSize,verbose=FALSE){
  s <- proc.time()
  analysedFile <- file(filename,"r")
  analysedLines <- NULL
  result <- NULL
  count <- 0
  header <- unlist(strsplit(readLines(analysedFile,n=1),"\t"))
  n.columns <- length(header)+1
  while(!((length(analysedLines) == 0) && (typeof(analysedLines) == "character"))){
    analysedLines<-readLines(con=analysedFile,n=sliceSize)
    if(!((length(analysedLines) == 0) && (typeof(analysedLines) == "character"))){
      aa <- t(matrix(unlist(strsplit(analysedLines,"\t")),n.columns,length(analysedLines)))
      analysedLines <- NULL
      gc()
      selected <- rbind(unlist(apply(aa,1,load.line.internal,pattern)))
      selected <- t(matrix(selected,n.columns,length(selected)/n.columns))
      result <- rbind(result,selected)
      count <- count + sliceSize
      if(verbose)cat("Processed: ",count,"lines\n")
    }
  }
  invisible(result)
}

read.populationHT.internal <- function(offspring,founders,map,founders_groups, populationType=c("riself", "f2", "bc", "risib"), threshold=0.01, sliceSize=5000, verbose=FALSE, debugMode=0){
  #**********FOUNDERS GROUPS*************
  if(missing(founders_groups)){ 
    stop("Specify founders_groups!\n")
  }else if(!(all(founders_groups %in% c(0,1)))){
    stop("founders_groups should contain only 0s and 1s.")
  }else if(all(founders_groups==1) || all(founders_groups==0)){
    stop("founders_groups should contain both 0s and 1s.")
  }
  #**********READING FOUNDERS PHENOTYPIC DATA*************
  filename <- paste(founders,"_phenotypes.txt",sep="")
  if(file.exists(filename)){
    population <- NULL
    if(verbose) cat("Found phenotypic file for founders:",filename,"and will store  it in population$founders$phenotypes\n")
    founders_phenotypes <- t.test.linebyline.internal(filename,founders_groups,threshold,sliceSize,verbose)
    population <- add.to.populationSub.internal(population,founders_phenotypes,"founders",populationType=populationType)
    population$founders$groups <- founders_groups
    doCleanUp.internal()
  }else{
    stop("No phenotype file for founders: ",filename," this file is essential in HT mode!\n")
  }
  #**********READING OFFSPRING PHENOTYPIC DATA*************
  filename <- paste(offspring,"_phenotypes.txt",sep="")
  if(file.exists(filename)){
    if(verbose) cat("Found phenotypic file for offspring:",filename,"and will store  it in population$offspring$phenotypes\n")
    offspring_phenotypes <- load.linebyline.internal(filename,rownames(population$founders$phenotypes),sliceSize,verbose)
    population <- add.to.populationSub.internal(population,offspring_phenotypes,"offspring_phenotypes",populationType=populationType)
    doCleanUp.internal()
  }else{
    stop("No phenotype file for offspring: ",filename," this file is essential!\n")
  }
  #**********READING CHILDREN GENOTYPIC DATA*************
  filename <- paste(offspring,"_genotypes.txt",sep="")
  if(file.exists(filename)){
    if(verbose) cat("Found genotypic file for offspring:",filename,"and will store  it in population$offspring$genotypes$real\n")
    offspring_genotypes <- load.linebyline.internal(filename,rownames(population$founders$phenotypes),sliceSize,verbose)
    population <- add.to.population(population, offspring_genotypes, "offspring$genotypes",populationType=populationType)
    doCleanUp.internal()
  }else{
    if(verbose)cat("No genotypic file for offspring:",filename,"genotypic data for offspring will be simulated\n")
  }
  #**********READING PHYSICAL MAP*************
  filename <- paste(map,"_physical.txt",sep="")
  if(file.exists(filename)){
    if(verbose) cat("Found genotypic file for offspring:",filename,"and will store  it in population$maps$physical\n")
    physical <- load.linebyline.internal(filename,rownames(population$founders$phenotypes),sliceSize,verbose)
    population <- add.to.population(population, physical, "maps$physical",populationType=populationType)
    doCleanUp.internal()
  }else{
    if(verbose)cat("No physical map file:",filename,".\n")
  }
  #**********READING GENETIC MAP*************
  filename <- paste(map,"_genetic.txt",sep="")
  if(file.exists(filename)){
    if(verbose) cat("Found genotypic file for offspring:",filename,"and will store  it in population$maps$genetic\n")
    genetic <- load.linebyline.internal(filename,rownames(population$founders$phenotypes),sliceSize,verbose)
    population <- add.to.population(population, genetic, "maps$genetic",populationType=populationType)
    doCleanUp.internal()
  }else{
    if(verbose)cat("No genetic map file:",filename,".\n")
  }
  invisible(population)
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
