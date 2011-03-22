#################################################################################
#
# parentalRoutine.R
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
# Contains: parentalRoutine 
# 				readParentalExpression, rankParentalExpression, filterParentalExpression
#
#################################################################################

#preprocessData
preprocessData <- function(expressionParental,groupLabels=c(0,0,1,1),verbose=FALSE,debugMode=0,...){
	s2<-proc.time()
	if(file.exists("RP.Rdata")){
		if(verbose) cat("File RP.Rdata already exists, reading it.\n")
		load("RP.Rdata")
	}else{
		rankParental <- invisible(RP(expressionParental,groupLabels,...))
		save(file="RP.Rdata",rankParental)
	}
	e2<-proc.time()
	if(verbose && debugMode==2)cat("Data preprocessing done in:",(e2-s2)[3],"seconds.\n")
	invisible(rankParental)
}


