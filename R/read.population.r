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
read.population <- function(offspring = "offspring", founders = "founders", map = "maps", founders_groups, populationType = c("riself", "f2", "bc", "risib"),
  read_mode = c("normal","HT"), verbose = FALSE, debugMode = 0,...){
  s <- proc.time()
  if(verbose && debugMode==1) cat("read.population starting.\n")
  population <- NULL
  populationType <- match.arg(populationType)
  read_mode <- match.arg(read_mode)
  if(read_mode=="normal"){
    population <- read.populationNormal.internal(offspring, founders, map, founders_groups, populationType, verbose, debugMode)
    population <- find.diff.expressed(population)
  }else{
    filename   <- paste(offspring,"_annotations.txt",sep="")
    population <- read.populationHT.internal(offspring, founders, map, founders_groups, populationType, verbose, debugMode,...)
    if(!file.exists(filename)){
      cat("No annotation file found in:",filename,"\n") #TODO: Warning / Error ?
      population <- find.diff.expressed(population)
    }
  }
  #**********FINALIZING FUNCTION*************
  e <- proc.time()
  if(verbose) cat("read.population finished after",(e-s)[3],"seconds.\n")
  doCleanUp.internal()
  invisible(population)
}

read.populationNormal.internal <- function(offspring,founders,map,founders_groups,populationType=c("riself", "f2", "bc", "risib"), verbose=FALSE,debugMode=0){
  #**********READING CHILDREN PHENOTYPIC DATA*************
  filename <- paste(offspring,"_phenotypes.txt",sep="")
  if(!file.exists(filename)) stop("No phenotype file for offspring: ",filename," this file is essentiall, you have to provide it\n")

  if(verbose) cat("Found phenotypic file for offspring:",filename,"and will store  it in population$offspring$phenotypes\n")
  offspring_phenotypes <- read.table(filename,sep="\t",header=TRUE)
  offspring_phenotypes <- as.matrix(offspring_phenotypes)
  population <- NULL
  population <- add.to.populationSub.internal(population,offspring_phenotypes,"offspring$phenotypes",populationType=populationType)
  doCleanUp.internal()
  
  #**********READING PARENTAL PHENOTYPIC DATA*************
  filename <- paste(founders,"_phenotypes.txt",sep="")
  if(missing(founders)){
    population <- simulateParentalPhentypes(population, population$offspring$phenotypes, populationType)
  }else if(file.exists(filename)){
    if(verbose) cat("Found phenotypic file for founders:",filename,"and will store it in population$founders$phenotypes\n")
    founders <- read.table(filename,sep="\t",header=TRUE)
    founders <- as.matrix(founders)
    population <- add.to.populationSub.internal(population, founders, "founders",populationType=populationType)
    #removing from founders probes that are not in children:
    population$founders$phenotypes <- mapMarkers.internal(population$founders$phenotypes,population$offspring$phenotypes, mapMode=1, verbose=verbose)
    population$founders$groups <- founders_groups
    doCleanUp.internal()
  }else{
    stop("No phenotype file for founders at location: ",filename," \n")
  }
  
  if(missing(population$founders$groups)) stop("Specify founders_groups!\n")

  #*****CHECK FOUNDERS GROUPS*************
  if(length(population$founders$groups) != ncol(population$founders$phenotypes)) stop("Length of founders_groups should be equal to the number of columns in founder phenotype file.")
  if((sum(population$founders$groups == 1) + sum(population$founders$groups == 0)) != length(founders_groups)) stop("founders_groups should contain only 0s and 1s.")

  class(population) <- c("population", populationType)
  #**********READING CHILDREN GENOTYPIC DATA*************
  filename <- paste(offspring,"_genotypes.txt",sep="")
  if(file.exists(filename)){
    if(verbose) cat("Found genotypic file for offspring:",filename,"and will store  it in population$offspring$genotypes$real\n")
    offspring_genotypes <- read.table(filename,sep="\t",header=TRUE)
    offspring_genotypes <- as.matrix(offspring_genotypes)
    population <- add.to.population(population, offspring_genotypes, "offspring$genotypes")
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
    population <- add.to.population(population, maps_genetic, "maps$genetic")
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
    population <- add.to.population(population, physical, "maps$physical")
    doCleanUp.internal()
  }else{
    if(verbose)cat("No physical map file:",filename,".\n")
  }
  invisible(population)
}

TwoClassTtest <- function(foundersline,foundersGroups,threshold){
  p1  <- as.numeric(foundersline[which(foundersGroups == 0) + 1])  #TODO: Why add 1, Why dont we get a out-of-bounds error ??
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
    cat("No phenotype file for founders at location: ", fileFoundersPheno," \n")
    #TODO: SIMULATE DATA WHEN NOT AVAILABLE !!!!!
  }
}

read.populationHT.internal <- function(offspring, founders, map, founders_groups, populationType=c("riself", "f2", "bc", "risib"), threshold=0.01,
  sliceSize = 5000, transformations = c("nothing","log","sqrt","reciprocal","probit","logit"), annots, verbose=FALSE, debugMode=0){
  #**********FOUNDERS GROUPS*************
  if(missing(founders_groups)){ 
    stop("Specify founders_groups!")
  }else if(!(all(founders_groups %in% c(0,1)))){
    stop("founders_groups should contain only 0s and 1s.")
  }else if(all(founders_groups==1) || all(founders_groups==0)){
    stop("founders_groups should contain both 0s and 1s.")
  }
  fileFoundersPheno  <- paste(founders,"_phenotypes.txt",sep="")
  fileOffspringPheno <- paste(offspring,"_phenotypes.txt",sep="")
  fileOffspringGeno  <- paste(offspring,"_genotypes.txt",sep="")
  fileAnnotations    <- paste(offspring,"_annotations.txt",sep="")
  fileMapPhys        <- paste(map,"_physical.txt",sep="")
  fileMapGen         <- paste(map,"_genetic.txt",sep="")
  population <- NULL
  if(!(file.exists(fileAnnotations))){

    #**********READING FOUNDERS PHENOTYPIC DATA*************
    if(verbose)cat("No annotation file",fileAnnotations,"\n")
    if(missing(founders)) stop("Founder phenotype data needs to be available:", fileFoundersPheno)
    readFoundersAndTtest(fileFoundersPheno, founders_groups, populationType, threshold, sliceSize, transformations, verbose)

    #**********READING OFFSPRING PHENOTYPIC DATA*************
    if(file.exists(fileOffspringPheno)){
      if(verbose) cat("Found phenotypic file for offspring:", fileOffspringPheno,"and will store  it in population$offspring$phenotypes\n")
      offspring_phenotypes <- loadByLine(fileOffspringPheno, rownames(population$founders$phenotypes), sliceSize, transformations, verbose)
      population <- add.to.populationSub.internal(population,offspring_phenotypes,"offspring$phenotypes",populationType=populationType)
      doCleanUp.internal()
    }else{
      stop("No phenotype file for offspring: ",fileOffspringPheno," this file is essential!\n")
    }
  }else{
    #**********READING FOUNDERS PHENOTYPIC DATA*************
    if(verbose) cat("Annotation file", fileAnnotations,"\n")
    population$flags  <- c(population$flags,"annots")
    population$annots <- read.table(fileAnnotations,sep="\t")
    readFoundersAndTtest(fileFoundersPheno, founders_groups, populationType, threshold, sliceSize, transformations, verbose)

    #**********READING OFFSPRING PHENOTYPIC DATA*************
    if(!file.exists(fileOffspringPheno)) stop("No phenotype file for offspring: ",fileOffspringPheno," this file is essential!\n")
    if(verbose) cat("Found phenotypic file for offspring:", fileOffspringPheno,"and will store  it in population$offspring$phenotypes\n")
    population$offspring$phenotypes <- fileOffspringPheno

  }
  class(population) <- c("population", populationType)
  #**********READING OFFSPRING GENOTYPIC DATA*************

  if(file.exists(fileOffspringGeno)){
    if(verbose) cat("Found genotypic file for offspring:",fileOffspringGeno,"and will store  it in population$offspring$genotypes$real\n")
    offspring_genotypes <- read.table(fileOffspringGeno,sep="\t", row.names=1, header=TRUE)
    population <- add.to.population(population, offspring_genotypes, "offspring$genotypes")
    doCleanUp.internal()
  }else{
    if(verbose) cat("No genotypic file for offspring:",fileOffspringGeno,"genotypic data for offspring will be simulated\n")
  }

  #**********READING PHYSICAL MAP*************
  if(file.exists(fileMapPhys)){
    if(verbose) cat("Found physical map:",fileMapPhys,"and will store  it in population$maps$physical\n")
    physical <- read.table(fileMapPhys,sep="\t", row.names=1, header=FALSE)
    population <- add.to.population(population, physical, "maps$physical")
    doCleanUp.internal()
  }else{
    if(verbose) cat("No physical map file:",fileMapPhys,".\n")
  }
  #**********READING GENETIC MAP*************
  if(file.exists(fileMapGen)){
    if(verbose) cat("Found genetic map:",fileMapGen,"and will store  it in population$maps$genetic\n")
    genetic <- read.table(fileMapGen,sep="\t",row.names=1,header=FALSE)
    population <- add.to.population(population, genetic, "maps$genetic")
    doCleanUp.internal()
  }else{
    if(verbose) cat("No genetic map file:",fileMapGen,".\n")
  }
  population$sliceSize <- sliceSize
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
