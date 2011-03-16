#################################################################################
#
# dataRoutines.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified March 2011
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
# Contains: readParentalExpression, readChildrenExpression, correctExpression
#           
#
#################################################################################


readParentalExpression <- function(filename,groups,treshold=0.05){
	#reading data in, user need to specify groups - sequence of 0s and 1s
	#this will be used in product rank function -> selected genes are loaded
}


readChildrenExpression <- function(filename,parentalOutput, correction=TRUE){
	#reading data in, corrects it and uses output from parental function
	# to select genes to further analysis
}


correctExpression <- function(){
	#using Danislav's batchcorrection probably, must think if it is wise to produce uber-massive
	#cross file...
}
