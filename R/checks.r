############################################################################################################
#
# checks.R
#
#
# Copyright (c) 2011, Konrad Zych
#
# 
# first written June 2011
# last modified July 2011
# last modified in version: 0.8.1
# in current version: active, in main workflow
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
# CONTAINS:
#	numericCheck.internal, genotypeCheck.internal, inRangeCheck.internal, inListCheck.internal, is.population
#
############################################################################################################

############################################################################################################
#									*** numericCheck.internal ***
#
# DESCRIPTION:
#	checking if given object is numeric or could be converted to numeric
# 
# PARAMETERS:
#	objectToBeChecked - object to be checked
#	allow.na - specifies whether object could cointain NAs
#
# OUTPUT:
#	boolean
#
############################################################################################################
numericCheck.internal <- function(objectToBeChecked, allow.na=FALSE){
	if(any(is.na(as.numeric(objectToBeChecked)))){
		### if object converted to numeric contains NAs it's either not-convertable or cointained NAs
		### previously, if we don't allow NAs - we just quit here, if we do, we have to chech if they
		### are intruduced (then object is not-convertable and we rerturm false) or they were in the
		### object previously (then we return true)
		if(!(allow.na)){
			return(FALSE)
		}else{
			if(sum(is.na(as.numeric(as.matrix(objectToBeChecked))))==sum(is.na((objectToBeChecked)))){
				return(TRUE)
			}else{
				return(FALSE)
			}
		}
	}
	return(TRUE)
}

############################################################################################################
#									*** genotypeCheck.internal ***
#
# DESCRIPTION:
#	checking if given object is containing only 0,1 and NAs
# 
# PARAMETERS:
# objectToBeChecked - object to be checked
# allow.na - specifies whether object could cointain NAs
#
# OUTPUT:
#	boolean
#
############################################################################################################
genotypeCheck.internal <- function(objectToBeChecked, allow.na=FALSE){
	converted <- as.numeric(as.matrix(objectToBeChecked))
	if(any(is.na(converted))){
		if(!(allow.na)){
			return(FALSE)
		}else{
			if(sum(is.na(converted))==sum(is.na(objectToBeChecked))){
				nrOfCorrect <- (sum(is.na(converted)) + sum(converted==1,na.rm=TRUE) + sum(converted==0,na.rm=TRUE))
				if(nrOfCorrect==length(converted)){
					return(TRUE)
				}else{
					return(FALSE)
				}
			}else{
				return(FALSE)
			}
		}
	}
	return(TRUE)
}

############################################################################################################
#									*** inRangeCheck.internal ***
#
# DESCRIPTION:
#	checking if given object is numeric and if it is from specified range
# 
# PARAMETERS:
# objectToBeChecked - object to be checked
# objectName - name of the object
# downLimit, upLimit - down and up limit to specify range
#
# OUTPUT:
#	none
#
############################################################################################################
inRangeCheck.internal <- function(objectToBeChecked, objectName, downLimit, upLimit){
	if(downLimit > upLimit){
		dl <- downLimit
		downLimit <- upLimit
		upLimit <- dl
	}
	if(downLimit==upLimit) warning("Specified range is just a single number!\n")
	if(!(numericCheck.internal(objectToBeChecked))) stop(objectName," is not numeric\n")
	if(objectToBeChecked<downLimit||objectToBeChecked>upLimit) stop(objectName," is: ",objectToBeChecked," but should be between ",downLimit," and ",upLimit,"\n")
}

############################################################################################################
#									*** inListCheck.internal ***
#
# DESCRIPTION:
#	checking if given object is numeric and if it is from specified range
# 
# PARAMETERS:
# objectToBeChecked - object to be checked
# listOfPossibleElements - 
#
# OUTPUT:
#	none
#
############################################################################################################
inListCheck.internal <- function(objectToBeChecked,objectName,listOfPossibleElements){
	if(length(listOfPossibleElements)==1) warning("Specified list contains just a single element!\n")
	if(any(!(objectToBeChecked%in%listOfPossibleElements))) stop("Selected wrong ",objectName," possibilities are: ",paste(listOfPossibleElements,sep=", "),"\n")
}

############################################################################################################
#									*** is.population ***
#
# DESCRIPTION:
#	checking if given object is of class population
# 
# PARAMETERS:
# objectToBeChecked - object to be checked
#
# OUTPUT:
#	none
#
############################################################################################################
is.population <- function(objectToBeChecked){
	if(is.null(objectToBeChecked)) stop("Supported object is empty!\n")
	if(is.null(class(objectToBeChecked)!="population")) stop("Supported object is not of a class population.\n")
	if(is.null(objectToBeChecked$offspring$phenotypes)) stop("No offspring phenotype data found, this is not a valid object of class population.\n")
	if(is.null(objectToBeChecked$founders$phenotypes)) stop("No founders phenotype data found, this is not a valid object of class population.\n")
	if(is.null(objectToBeChecked$founders$groups)) stop("No information about founders groups found, this is not a valid object of class population.\n")
}

############################################################################################################
#									*** crossContainsMap.internal  ***
#
# DESCRIPTION:
# 	checking if crss object contains specified type of map
# 
# PARAMETERS:
# 	cross - object of class cross, containing physical or genetic map
#	map - which map should be used for comparison:
#			- genetic - genetic map from cross$maps$genetic
#			- physical - physical map from cross$maps$physical
# 
# OUTPUT:
#	none
#
############################################################################################################
crossContainsMap.internal <- function(cross,map=c("genetic","physical")){
	if(map=="genetic"){
		if(is.null(cross$maps$genetic)){
			stop("Selected map=genetic, but there is no genetic map in cross$maps$genetic")
		}
	}else if(map=="physical"){
		if(is.null(cross$maps$physical)){
			stop("Selected map=physical, but there is no physical map in cross$maps$physical")
		}
	}else{
		stop("Parameter map should be either genetic or physical")
	}
}

############################################################################################################
#									*** defaultCheck.internal  ***
#
# DESCRIPTION:
# 	making sure that default parameter is used, when parameter is speicified by =c("","")
# 
# PARAMETERS:
# 	parameterToBeChecked
#	maxLength
# 
# OUTPUT:
#	default parameter from list of possible
#
############################################################################################################
defaultCheck.internal <- function(parameterToBeChecked,nameOfParameter,maxLength){
	if(length(parameterToBeChecked)==maxLength){
		parameterToBeChecked <- parameterToBeChecked[1]
	}else if(length(parameterToBeChecked)!=1){
		stop("wrong parameter ",nameOfParameter," length, choose one out of possible\n")
	}
}
