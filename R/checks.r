#
# checks.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Apr, 2013
# first written Dec, 2011
# Contains: numericCheck.internal, genotypeCheck.internal, inRangeCheck.internal 
#           inListCheck.internal, check.population, crossContainsMap.internal 
#           defaultCheck.internal
#

# numericCheck.internal
#
# DESCRIPTION:
#  Checking if given object is numeric or could be converted to numeric
# OUTPUT:
#  boolean
#
numericCheck.internal <- function(objectToBeChecked, allow.na=FALSE){
  converted <- as.numeric(as.matrix(objectToBeChecked))
  if(any(is.na(converted))){
    if(!(allow.na)) return(FALSE)
    if(sum(is.na(converted)) == sum(is.na((objectToBeChecked)))) return(TRUE)
    return(FALSE)
  }
  return(TRUE)
}

# genotypeCheck.internal
#
# DESCRIPTION:
#  checking if given object is containing only 0,1 and NAs
# OUTPUT:
#  boolean
#
genotypeCheck.internal <- function(objectToBeChecked, genotypes, allow.na=FALSE){
  converted <- as.numeric(as.matrix(objectToBeChecked))
  if(any(is.na(converted))){
    if(!(allow.na)) return(FALSE)
    if(sum(is.na(converted))==sum(is.na(objectToBeChecked))){
      nrOfCorrect <- sum(is.na(converted))
      for(genotype in genotypes){ nrOfCorrect <- nrOfCorrect + sum(converted==genotype,na.rm=TRUE) }
      if(nrOfCorrect==length(converted)) return(TRUE)
      return(FALSE)
    }
    return(FALSE)
  }
  return(TRUE)
}

############################################################################################################
#                                       *** inRangeCheck.internal ***
#
# DESCRIPTION:
#  checking if given object is numeric and if it is from specified range
# OUTPUT:
#  none
############################################################################################################
inRangeCheck.internal <- function(objectToBeChecked, objectName, downLimit, upLimit){
  if(downLimit > upLimit){
    dl <- downLimit
    downLimit <- upLimit
    upLimit <- dl
  }
  if(!(numericCheck.internal(objectToBeChecked))) stop(objectName," is not numeric\n")
  if(objectToBeChecked<downLimit||objectToBeChecked>upLimit) stop(objectName," is: ",objectToBeChecked," but should be between ",downLimit," and ",upLimit,"\n")
}

############################################################################################################
#                  *                     ** inListCheck.internal ***
#
# DESCRIPTION:
#  checking if given object is numeric and if it is from specified range
# OUTPUT:
#  none
############################################################################################################
inListCheck.internal <- function(objectToBeChecked,objectName,listOfPossibleElements){
  if(length(listOfPossibleElements)==1) warning("Specified list contains just a single element!\n")
  if(any(!(objectToBeChecked%in%listOfPossibleElements))) stop("Selected wrong ",objectName," possibilities are: ",paste(listOfPossibleElements,sep=", "),"\n")
}

############################################################################################################
#                                          *** check.population ***
#
# DESCRIPTION:
#  checking if given object is a correct population object
# OUTPUT:
#  none
############################################################################################################
check.population <- function(x,verbose=FALSE){
  #if(length(class(x))!=2) stop("Incorrect class of the object.\n")
  if(class(x)[1]!="population") stop("Object is not of a class population.\n")
  if(!(class(x)[2] %in% c("riself", "f2", "bc", "risib"))) stop("Type of the population: ",class(x)[2]," not recognized.\n")
  if(is.null(x$founders$phenotypes) && !("noParents" %in% x$flags) && !("annots" %in% x$flags)) stop("No offspring phenotype data found, this is not a valid object of class population.\n")
  if(is.null(x$offspring$phenotypes)) stop("No offspring phenotype data found, this is not a valid object of class population.\n")
  if(is.null(x$founders$groups) && !("noParents" %in% x$flags) && !("annots" %in% x$flags)) stop("No offspring phenotype data found, this is not a valid object of class population.\n")
}

############################################################################################################
#                                      *** checkParameters.internal  ***
#
# DESCRIPTION:
#   making sure that default parameter is used, when parameter is speicified by =c("","")
# OUTPUT:
#  default parameter from list of possible
#
############################################################################################################
checkParameters.internal <- function(parameterToBeChecked,possibleValues,nameOfParameter=""){
  if(length(parameterToBeChecked)==length(possibleValues)){
    if(any(!(parameterToBeChecked%in%possibleValues))) stop(nameOfParameter," parameter is incorrect, possible values: ",paste(possibleValues,sep="\t"),"\n")
    return(parameterToBeChecked[1])
  }else if(length(parameterToBeChecked)==1){
    if(!(parameterToBeChecked%in%possibleValues)) stop(nameOfParameter," parameter is incorrect, possible values: ",paste(possibleValues,sep="\t"),"\n")
    return(parameterToBeChecked)
  }else{
    stop(nameOfParameter," parameter is incorrect, possible values: ",paste(possibleValues,sep="\t"),"\n")
  }
}

############################################################################################################
#                  *** defaultCheck.internal  ***
#
# DESCRIPTION:
#   making sure that default parameter is used, when parameter is speicified by =c("","")
# OUTPUT:
#  default parameter from list of possible
############################################################################################################
defaultCheck.internal <- function(parameterToBeChecked,nameOfParameter,maxLength,defVal){
  if(length(parameterToBeChecked) == maxLength){
    invisible(defVal)
  }else if(length(parameterToBeChecked) != 1){
    stop("wrong parameter ",nameOfParameter," length, choose one out of possible\n")
  }else{
    invisible(parameterToBeChecked)
  }
}

