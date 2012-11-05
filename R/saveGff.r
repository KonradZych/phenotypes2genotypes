#
# saveGff.r
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Nov, 2012
# first written Nov, 2012
# Contains: saveGff
#

#  saveGff
#
# DESCRIPTION:
#  Saves a gff file containing physical locations for markers for use in genome viewers.
# PARAMETERS:
#   - population -  An object of class population.
#   - file - an output file
#   - verbose - Be verbose
#   - debugMode - 1: Print our checks, 2: print additional time information
# OUTPUT:
#  NONE
#

saveGff.internal <- function(gffFile="population.gff", markers, positions){
  cat("##gff-version 3\n",file=gffFile,append=FALSE)
  for(marker in rownames(markers)){
    cat(marker,"\t.\tmarker\t",positions[marker,],"\t",positions[marker,],"\t.\t+\t.\tID=",marker,"\n",file=gffFile,append=TRUE)
  }
}