#
# genotypesToCross.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Mar, 2011
# Contains: genotypesToCross.internal, writePhenotypes.internal, 
#           writeGenotypes.internal, cleanNames.internal
#

# genotypesToCross.internal
#
# DESCRIPTION:
#  Produces from genotypic matrix file containing object of type cross, reads it into R a returns
# PARAMETERS:
#   - population - population type object, must contain founders phenotypic data.
#   - use - save "real" gentypes, "simulated" genotypes otr simulated genotypes ordered using "map" from gff file
#   - outputFile - file where object of type cross is being saved
#   - verbose - Be verbose
#   - debugMode - 1: Print our checks, 2: print additional time information
# OUTPUT:
#  An object of class cross
#
genotypesToCross.internal <- function(population, genotype=c("simulated","real"), orderUsing=c("none","map_genetic","map_physical"), outputFile="mycross.csv", verbose=FALSE, debugMode=0){
  #checks
  if(missing(population)) stop("No population object provided.\n") 
  check.population(population)
  orderUsing <- checkParameters.internal(orderUsing,c("none","map_genetic","map_physical"),"orderUsing")
  comparisonMethod <- checkParameters.internal(genotype,c("simulated","real"),"genotype")
  if(orderUsing=="map_physical"&&is.null(population$maps$physical)) stop("There is no map in population$maps$physical\n")
  if(orderUsing=="map_genetic"&&is.null(population$maps$genetic)) stop("There is no map in population$maps$genetic\n")
  if(verbose && debugMode==1) cat("genotypesToCross starting without errors in checkpotins.\n")
  s <- proc.time()
  
  #WRITING PHENOTYPIC DATA TO FILE
  population<-writePhenotypes.internal(population, genotype, outputFile, verbose, debugMode)
  
  #WRITING GENOTYPIC DATA TO FILE
  if(genotype=="real"){
    if(is.null(population$offspring$genotypes$real)){
      stop("Use = real chosen, but there is no real genotypic data in population$offspring$genotypes$real\n")
    }else{
      genoL <- length(table(population$offspring$genotypes$real))
      if(orderUsing=="none"){
        writeGenotypes.internal(population$offspring$genotypes$real, chr=1, outputFile=outputFile, verbose=verbose, debugMode=debugMode)
        genotypes <- names(table(population$offspring$genotypes$real))
      }else if(orderUsing=="map_physical"){
        population$maps$physical <- mapMarkers.internal(population$maps$physical,population$offspring$genotypes$real, mapMode=1, verbose=verbose)
        if(is.null(population$maps$physical)) stop("No physical map provided in population$maps$physical\n")
        writeGenotypes.internal(population$offspring$genotypes$real, chr=population$maps$physical[rownames(population$offspring$genotypes$real),1], positions=population$maps$physical[rownames(population$offspring$genotypes$real),2], outputFile=outputFile, verbose=verbose, debugMode=debugMode)
        genotypes <- names(table(population$offspring$genotypes$real))
      }else if(orderUsing=="map_genetic"){
        population$maps$genetic <- mapMarkers.internal(population$maps$genetic,population$offspring$genotypes$real, mapMode=1, verbose=verbose)
        if(is.null(population$maps$genetic)) stop("No genetic map provided in population$maps$genetic\n")
        writeGenotypes.internal(population$offspring$genotypes$real, chr=population$maps$genetic[rownames(population$offspring$genotypes$real),1], positions=population$maps$genetic[rownames(population$offspring$genotypes$real),2], outputFile=outputFile, verbose=verbose, debugMode=debugMode)
        genotypes <- names(table(population$offspring$genotypes$real))
      }
    
    }
  }else if(genotype=="simulated"){
    if(is.null(population$offspring$genotypes$simulated)){
      stop("Use = simulated chosen, but there is no simulated genotypic data in population$offspring$genotypes$simulated\n")
    }else{
      genoL <- length(table(population$offspring$genotypes$simulated))
      if(orderUsing=="none"){
        writeGenotypes.internal(population$offspring$genotypes$simulated, chr=1, outputFile=outputFile, verbose=verbose, debugMode=debugMode)
        genotypes <- names(table(population$offspring$genotypes$simulated))
      }else if(orderUsing=="map_physical"){
        population$maps$physical <- mapMarkers.internal(population$maps$physical,population$offspring$genotypes$simulated, mapMode=1, verbose=verbose)
        if(is.null(population$maps$physical)) stop("orderUsing = map_physical chosen, but no physical map provided in population$maps$physical\n")
        writeGenotypes.internal(population$offspring$genotypes$simulated, chr=population$maps$physical[rownames(population$offspring$genotypes$simulated),1], positions=population$maps$physical[rownames(population$offspring$genotypes$simulated),2], outputFile=outputFile, verbose=verbose, debugMode=debugMode)
        genotypes <- names(table(population$offspring$genotypes$simulated))
      }else if(orderUsing=="map_genetic"){
        population$maps$genetic <- mapMarkers.internal(population$maps$genetic,population$offspring$genotypes$simulated, mapMode=1, verbose=verbose)
        if(is.null(population$maps$genetic)) stop("orderUsing = map_physical chosen, but no genetic map provided in population$maps$genetic\n")
        writeGenotypes.internal(population$offspring$genotypes$simulated, chr=population$maps$genetic[rownames(population$offspring$genotypes$simulated),1], positions=population$maps$genetic[rownames(population$offspring$genotypes$simulated),2], outputFile=outputFile, verbose=verbose, debugMode=debugMode)
        genotypes <- names(table(population$offspring$genotypes$simulated))
      }
    }
  }  

  #READING CROSS FILE INTO R
  cross <- invisible(read.cross("csvr",file=outputFile, genotypes=genotypes))
  cross <- convertType.internal(cross,class(population)[2])
  e <- proc.time()
  if(verbose) cat("genotypesToCross done in",(e-s)[3],"seconds.\n")
  invisible(cross)
}

convertType.internal <- function(cross,populationType){
	if(populationType == "riself"){
		cross <- convert2riself(cross)
	}else if(populationType == "risib"){
		cross <- convert2risib(cross)
	}else{
		class(cross)[1] <- populationType
	}
	return(cross)
}

############################################################################################################
#                  *** writePhenotypes.internal ***
#
# DESCRIPTION:
#  sub function of genotypesToCross - writes phenotypes to file
# 
# PARAMETERS:
#   population - Ril type object, must contain founders phenotypic data.
#   use - save "real" gentypes, "simulated" genotypes otr simulated genotypes ordered using "map" from gff file
#   outputFile - file where object of type cross is being saved
#   verbose - Be verbose
#   debugMode - 1: Print our checks, 2: print additional time information
#
# OUTPUT:
#  none
#
############################################################################################################
writePhenotypes.internal <- function(population, genotype, outputFile, verbose=FALSE, debugMode=0){
  sl <- proc.time()
  if(verbose && debugMode==1) cat("writePhenotypes starting.\n")
  if(genotype=="real"){
    if(is.null(population$offspring$genotypes$real)){
      stop("genotype = real chosen, but there is no real genotypic data in population$offspring$genotypes$read\n")
    }else{
      nr_b <- ncol(population$offspring$phenotypes)
      population$offspring$phenotypes <- mapMarkers.internal(population$offspring$phenotypes,population$offspring$genotypes$real,mapMode=2)
      population$offspring$genotypes$real <- mapMarkers.internal(population$offspring$genotypes$real,population$offspring$phenotypes,mapMode=2)
      nr_a <- ncol(population$offspring$phenotypes)
      if(verbose)cat(nr_b-nr_a,"individuals out of",nr_b,"were removed due to mismatch\n")
    }
  }else if(genotype=="simulated"){
    if(is.null(population$offspring$genotypes$simulated)){
      stop("genotype = simulated or map chosen, but there is no simulated genotypic data in population$offspring$genotypes$simulated\n")
    }else{
      nr_b <- ncol(population$offspring$phenotypes)
      population$offspring$phenotypes <- mapMarkers.internal(population$offspring$phenotypes,population$offspring$genotypes$simulated,mapMode=2)
      population$offspring$genotypes$simulated <- mapMarkers.internal(population$offspring$genotypes$simulated,population$offspring$phenotypes,mapMode=2)
      nr_a <- ncol(population$offspring$phenotypes)
      if(verbose)cat(nr_b-nr_a,"individuals out of",nr_b,"were removed due to mismatch\n")
    }
  }
  population$offspring$phenotypes<- cleanNames.internal(population$offspring$phenotypes)
  write.table(cbind("","",population$offspring$phenotypes),file=outputFile,sep=",",quote=FALSE,col.names=FALSE)
  el <- proc.time()
  if(verbose && debugMode==2)cat("Writing phenotypes done in:",(el-sl)[3],"seconds.\n")
  invisible(population)
}

############################################################################################################
#                  *** writeGenotypes.internal ***
#
# DESCRIPTION:
#  sub function of genotypesToCross and writeGenotypesMap - writes genotypes (one chromosome at the time) 
#  to file
# 
# PARAMETERS:
#   genotypeMatrix - matrix of genotypic data, rows - markers, cols - individuals
#   chr - chromosome currently being written
#   outputFile - file where object of type cross is being saved
#   verbose - Be verbose
#   debugMode - 1: Print our checks, 2: print additional time information
#
# OUTPUT:
#  none
#
############################################################################################################
writeGenotypes.internal <- function(genotypeMatrix,chr=1,positions=NULL,outputFile,verbose=FALSE,debugMode=0){
  sl <- proc.time()
  if(verbose && debugMode==1) cat("writeGenotypes starting.\n")
  if(is.null(positions)) positions <- 1:nrow(genotypeMatrix)
  else if(length(positions)!=length(1:nrow(genotypeMatrix))) stop("Posistions object is not correct, check help files.\n")
  if(verbose && debugMode==1) cat("writeGenotypes starting.\n")
  genotypeMatrix <- cleanNames.internal(genotypeMatrix)
  write.table(cbind(rownames(genotypeMatrix),chr,positions,genotypeMatrix),file=outputFile,sep=",",quote=FALSE,
    col.names=FALSE,append=TRUE,row.names=FALSE)
  el <- proc.time()
  if(verbose && debugMode==2) cat("Writing genotypes done in:",(el-sl)[3],"seconds.\n")
}

############################################################################################################
#                  *** cleanNames.internal ***
#
# DESCRIPTION:
#  changing names that will crush read.cross
# 
# PARAMETERS:
#   matrixToBeCleaned - matrix of any data type
#
# OUTPUT:
#  matrix of any data type
#
############################################################################################################
cleanNames.internal <- function(matrixToBeCleaned){
  for(i in 1:nrow(matrixToBeCleaned)){
    old <- rownames(matrixToBeCleaned)[i]
    new <- gsub(",","_",rownames(matrixToBeCleaned)[i])
    if(old != new){
      rownames(matrixToBeCleaned)[i] <- new
      cat("WARNING: marker name switched from:",old,"to",new,"because it contained ','!\n")
    }
  }
  invisible(matrixToBeCleaned)
}
