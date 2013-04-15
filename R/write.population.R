#
# write.population.r
#
# Copyright (c) 2010-2013 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Apr, 2013
# first written Apr, 2013
# Contains: write.population, writeSingleFile
#

#  write.population
#
# DESCRIPTION:
#  Writes geno/phenotypic files from R environment special object to disk.
# PARAMETERS:
#   - population - An object of class population .
#   - offspring - Core used to specify names of children phenotypic ("offspring_phenotypes.txt") and genotypic ("offspring_genotypes.txt") files.
#   - founders - Core used to specify names of founders phenotypic ("founders_phenotypes.txt") file.
#   - map - Core used to specify names of genetic ("map_genetic.txt") and physical ("map_physical.txt") map files.
#   - verbose - Be verbose
#   - debugMode - 1: Print our checks, 2: print additional time information
# OUTPUT:
#   None
#
write.population <- function(population, offspring = "offspring", founders = "founders", map = "map", verbose = FALSE, debugMode = 0){

  ### checks
  check.population(population)

  ### file names
  fileFoundersPheno  <- paste(founders,"_phenotypes.txt",sep="")
  fileOffspringPheno <- paste(offspring,"_phenotypes.txt",sep="")
  fileOffspringGeno  <- paste(offspring,"_genotypes.txt",sep="")
  fileAnnotations    <- paste(offspring,"_annotations.txt",sep="")
  fileMapPhys        <- paste(map,"_physical.txt",sep="")
  fileMapGen         <- paste(map,"_genetic.txt",sep="")

  ### initializing  
  s <- proc.time()
  
  ### offspring phenotypic file
  if(!is.null(dim(population$offspring$phenotypes))){
    writeSingleFile(population$offspring$phenotypes, "offspring phenotypes", fileOffspringPheno, verbose=verbose)
  }else{
    stop("Phenotype data for offspring are already stored in:",population$offspring$phenotypes,". The function will not overwrite that file.\n")
  }

  ### founders phenotypic file
  writeSingleFile(population$founders$phenotypes, "founder phenotypes", fileFoundersPheno, verbose=verbose)
  
  ### annotations file
  writeSingleFile(population$annots, "annotation", fileAnnotations, verbose=verbose)

  ### offspring genotypic file
  writeSingleFile(population$offspring$genotypes, "offspring genotypes", fileOffspringGeno, verbose=verbose)
  
  ### physical map
  writeSingleFile(population$maps$physical, "physical map", fileMapPhys, verbose=verbose, col.names=FALSE)

  ### genetic map
  writeSingleFile(population$maps$genetic, "genetic map", fileMapGen, verbose=verbose, col.names=FALSE)

  #**********FINALIZING FUNCTION*************
  e <- proc.time()
  if(verbose && debugMode==2) cat("read.population finished after",(e-s)[3],"seconds.\n")
}

#  writeSingleFile
#
# DESCRIPTION:
#  Writes single geno/phenotypic file.
# PARAMETERS:
#   - dataMatrix - an object of class population
#   - filename - name of the file that will be processed
#   - verbose - be verbose
#   - ... - passed to write.table
# OUTPUT:
#   None
#
writeSingleFile   <- function(dataMatrix, dataType, filename, verbose=FALSE, ...){
  if(file.exists(filename))  stop("File: ",filename," already exists!\n")
  if(!is.null(dim(dataMatrix))){
    write.table(dataMatrix,file=filename,sep="\t",quote=FALSE,...)
    if(verbose) cat(dataType,"saved in file:",filename,"\n")
  }else{
    if(verbose) cat("no",dataType,"found\n")
  }
}
