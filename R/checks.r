############################################################################################################
#
# checks.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written June 2011
# last modified June 2011
# last modified in version: 0.7.1
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
# Contains: 
#
############################################################################################################

############################################################################################################
#numericCheck.internal - checking if given object is numeric or could be converted to numeric
# 
# objectToBeChecked - object to be checked
# allow.na - specifies whether object could cointain NAs
# value - boolean
#
############################################################################################################
numericCheck.internal <- function(objectToBeChecked, allow.na=FALSE){
	if(any(is.na(as.numeric(objectToBeChecked)))){
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
#inRangeCheck.internal - checking if given object is numeric and if it is from specified range
# 
# objectToBeChecked - object to be checked
# objectName - name of the object
# downLimit, upLimit - down and up limit to specify range
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
#inRangeCheck.internal - checking if given object is numeric and if it is from specified range
# 
# objectToBeChecked - object to be checked
# objectName - name of the object
# downLimit, upLimit - down and up limit to specify range
#
############################################################################################################
inListCheck.internal <- function(objectToBeChecked,objectName,listOfPossibleElements){
	if(length(listOfPossibleElements)==1) warning("Specified list contains just a single element!\n")
	if(!(any(listOfPossibleElements==objectToBeChecked))) stop("Selected wrong ",objectName," possibilities are: ",listOfPossibleElements,"\n")
}
