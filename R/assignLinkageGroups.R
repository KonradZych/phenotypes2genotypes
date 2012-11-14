#
# assignLinkageGroups.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Oct, 2011
# Contains: assignLinkageGroups, regorganizeMarkersWithin
#

# assignLinkageGroups
#
# DESCRIPTION:
#  Assign linkage groups based on a user supplied known number of chromosomes
#  we can use the genetic map of a cross, or the rf matrix
# OUTPUT:
#  an object of class cross
#
assignLinkageGroups <- function(cross, n.chr, use=c("geno","rf"), ...){
  use <- match.arg(use)
  geno <- t(pull.geno(cross))
  inplaceOfNA <- min(geno,na.rm=T)-1
  geno[which(is.na(geno))] <- inplaceOfNA
  if(use=="geno") clustering <- kmeans(geno, n.chr, nstart=100, ...)
  if(use=="rf"){
    cross <- cleanRfs.internal(cross)
    clustering <- kmeans(lowerTrng.internal(est.rf(cross)$rf), n.chr, nstart=100, ...)
  }
  reorganizeMarkersWithin(cross, clustering$cluster)
}

cleanRfs.internal <- function(cross){
  dataRf <- est.rf(cross)$rf
  dataRf <- lowerTrng.internal(dataRf)
  minis <- apply(dataRf,1,min)
  for(i in 1:nrow(dataRf)){
    if(minis[i]>0.5){
      cat("dropping",rownames(dataRf)[i],"min:",minis[i])
      cross <- drop.markers(cross, rownames(dataRf)[i])
    }
  }
  invisible(cross)
}

lowerTrng.internal <- function(dataRf){
  for(i in 1:nrow(dataRf)){
    for(j in 1:i){
      dataRf[j,i] <- dataRf[i,j]
    }
  }
  invisible(dataRf)
}

# reorganizeMarkersWithin
#
# DESCRIPTION:
#  Function to quickly rearrange all the markers in a cross based on any ordering vector
# OUTPUT:
#  an object of class cross
#
reorganizeMarkersWithin <- function(cross, ordering){
  cross <- clean(cross)
  n.markers <- sum(nmar(cross))
  chrtype <- rep(sapply(cross$geno, class), n.markers)
  crosstype <- class(cross)[1]
  g <- pull.geno(cross)
  newChromosomes <- sort(unique(ordering))
  n.newChromosomes <- length(newChromosomes)
  cross$geno <- vector("list", n.newChromosomes)
  names(cross$geno) <- newChromosomes
  for (i in 1:n.newChromosomes) {
      selectedMarkers <- which(ordering == newChromosomes[i])
      cross$geno[[i]]$data <- g[, selectedMarkers, drop = FALSE]
      cross$geno[[i]]$map <- seq(0, by = 10, length = length(selectedMarkers))
      if (crosstype == "4way") {
          cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map, cross$geno[[i]]$map)
          colnames(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
      }
      else names(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
      thechrtype <- unique(chrtype[which(ordering == newChromosomes[i])])
      if (length(thechrtype) > 1) 
          warning("Problem with linkage group ", i, ": A or X?\n", 
            paste(thechrtype, collapse = " "))
      else class(cross$geno[[i]]) <- thechrtype
  }
  cross <- est.rf(cross)
  return(cross)
}
