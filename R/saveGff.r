#
# saveGff.r
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified Nov, 2012
# first written Nov, 2012
# Contains: saveGff, findrecombinations
#

#  saveGff
#
# DESCRIPTION:
#  Saves a gff file containing physical locations for markers for use in genome viewers.
# PARAMETERS:
#   - population -  An object of class cross.
#   - population -  An object of class population.
#   - gffFileCore - an output file
#   - verbose - Be verbose
# OUTPUT:
#  NONE
#
saveGff <- function(cross, population, gffFileCore="population", verbose=FALSE){
  if(missing(population)) stop("Please provide a population object\n")
  populationType <- class(population)[2]
  if(is.null(population$maps$physical)) stop("No physical map in the population object!\n")
  check.population(population)
  markersCross <- markernames(cross)
  markers <- markersCross[which(markersCross %in% rownames(population$maps$physical))]
  markers <- c(markers,cross$redundant)
  if(is.integer0(markers) || is.null(markers)) stop("No physical locations for any of the markers in the cross object!\n")
  crossNr <- drop.dupmarkers(cross,verbose=FALSE)
  markersCrossNr <- markernames(crossNr)
  markersNr <- markersCrossNr[which(markersCrossNr %in% rownames(population$maps$physical))]
  if(is.integer0(markersNr) || is.null(markersNr)) stop("No physical locations for any of the non redundant markers in the cross object!\n")
  recombinations <- findrecombinations(markers)
  recombinationsNr <- findrecombinations(markersNr)
  filename1 <- paste(gffFileCore,"_markers_nonRedundant.gff",sep="")
  filename2 <- paste(gffFileCore,"_recombinations_nonRedundant.gff",sep="")
  filename3 <- paste(gffFileCore,"_markers_all.gff",sep="")
  filename4 <- paste(gffFileCore,"_recombinations_all.gff",sep="")
  cat("##gff-version 3\n",file=filename1)
  cat("##gff-version 3\n",file=filename2)
  cat("##gff-version 3\n",file=filename3)
  cat("##gff-version 3\n",file=filename4)
  for(marker in markers){
      if(!(marker %in% markersNr)){
        cat("Chr",population$maps$physical[marker,1],"\t.\tmarker-red\t",population$maps$physical[marker,2],"\t",population$maps$physical[marker,3],"\t100\t+\t.\tID=",marker,"\n",file=filename3,append=TRUE,sep='')
      }else{
        cat("Chr",population$maps$physical[marker,1],"\t.\tmarker-red\t",population$maps$physical[marker,2],"\t",population$maps$physical[marker,3],"\t100\t+\t.\tID=",marker,"\n",file=filename1,append=TRUE,sep='')
        cat("Chr",population$maps$physical[marker,1],"\t.\tmarker-nonred\t",population$maps$physical[marker,2],"\t",population$maps$physical[marker,3],"\t100\t+\t.\tID=",marker,"\n",file=filename3,append=TRUE,sep='')
      }
  }
  if(verbose) cat("Saved ",length(markersNr),"non-redundant and",length(markers)-length(markersNr),"redundant markers. In total:",length(markers),"markers.\n")
}

findrecombinations <- function(markers){
  return(markers)
}
