#
# markersCorPlot.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Mar, 2011
# Contains: markersCorPlot, ascendingMaptoJigSawMap, getChrOffsets.internal
#           getMarkerOffsets, getMarkerOffsetsFromMap, chromCorMatrix
#           getPopulationOffsets.internal 
#

# markersCorPlot
#
# DESCRIPTION:
#  function to create new map and save it in cross object
# PARAMETERS:
#   - population - object of class population
#   - orde - object of class population
#   - n.chr - expected number of linkage groups
#   - use - expected number of linkage groups
#   - verbose - be verbose
# OUTPUT:
#  An object of class cross
#
qtl.comparison.plot <- function(cross, cross2, map.physical1, map.physical2, chr, ...){
  if(missing(cross))   stop("Provide two objects of class cross!\n")
  if(missing(cross2))  stop("Provide two objects of class cross!\n")
  if(!(any(class(cross)=="cross")))   stop("Provide two objects of class cross!\n")
  if(!(any(class(cross2)=="cross")))  stop("Provide two objects of class cross!\n")
  res1 <- scanone(cross,...)
  res2 <- scanone(cross2,...)
  if((missing(chr))){
    chr <- 1
  }
  if(length(chr)>1){
    chr <- chr[1]
    warning("Function can plot only one chromosome at the time, plotting first of selected chromosomes!")
  }
  if(!(chr%in%1:min(nchr(cross),nchr(cross2)))){
    stop("Chromosome asked is not present in the cross object\n")
  }
  if(!(missing(map.physical1))){
    #op <-par(mfrow=c(floor(sqrt(length(chr))),round(sqrt(length(chr)))))
    for(n in chr){
      max.lod1 <- max(res1[which(res1[,1]==n),3])
      max.lod2 <- max(res2[which(res2[,1]==n),3])
      max.lod <- max(max.lod1,max.lod2)
      plot(res1,res2,col=c("orange","black"),chr=n,main="Figure 2 - QTL profiles comparison.")
      values = NULL
      for(x in 1:(min(nchr(cross),nchr(cross2)))){
        values = c(values,max(map.physical1[which(map.physical1[,1]==x),2]))
      }
      scaler <- sum(values)/max.lod
      valuesS = 0
      chrPoints = NULL
      map.physical_ <- map.physical1
      map.physical_2 <- map.physical2
      for(x in 1:(min(nchr(cross),nchr(cross2)))){
        map.physical_[which(map.physical_[,1]==x),2] = map.physical1[which(map.physical1[,1]==x),2]+valuesS
        map.physical_2[which(map.physical_2[,1]==x),2] = map.physical2[which(map.physical2[,1]==x),2]+valuesS
        oldV = valuesS
        valuesS = valuesS+values[x]
        chrPoints = c(chrPoints,(oldV+valuesS)/2)
        abline(h=valuesS/scaler,lty=2)
      }
      chrPoints = c(chrPoints,(sum(values)+valuesS)/2)
      usefulMnames <- which(rownames(map.physical_)%in%markernames(cross))
      usefulMnames2 <- which(rownames(map.physical_2)%in%markernames(cross))
      usefulM <- map.physical_[usefulMnames,]
      usefulM2 <- map.physical_2[usefulMnames2,]
      usefulmonNewMAP<-which(as.character(names(cross$geno[[n]]$map))%in%as.character(rownames(usefulM)))
      usefulmonNewMAP2<-which(as.character(names(cross$geno[[n]]$map))%in%as.character(rownames(usefulM2)))
      usefulMChr <- usefulM[names(cross$geno[[n]]$map)[usefulmonNewMAP],]
      usefulMChr2 <- usefulM2[names(cross$geno[[n]]$map)[usefulmonNewMAP2],]
      points(cross$geno[[n]]$map[usefulmonNewMAP],(usefulMChr[,2]/scaler),pch=20,col="orange")
      points(cross$geno[[n]]$map[usefulmonNewMAP2],(usefulMChr2[,2]/scaler),pch=20)
      axis(4,chrPoints[1:(length(chrPoints)-1)]/scaler,paste("Chr",1:(min(nchr(cross),nchr(cross2))),sep=" "))
    }
  }
}