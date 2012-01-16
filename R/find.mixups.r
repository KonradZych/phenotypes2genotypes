############################################################################################################
#
# find.mixups.R
#
# Copyright (c) 2012, Konrad Zych
#
# Modified by Danny Arends
# 
# first written January 2012
# last modified January 2012
# last modified in version: 1.0.0
# in current version: active, in main workflow
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
#
#     This program is distributed in the hope that it will be useful
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
#
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: cross.saturate, rearrangeMarkers,
#             bestCorelated.internal, map2mapCorrelationMatrix.internal, map2mapImage
#
############################################################################################################

###########################################################################################################
#                                    *** find.mixups ***
#
# DESCRIPTION:
# 	saturate existing genetic map adding markers derived from gene expression
# OUTPUT:
#	object of class cross
############################################################################################################
find.mixups <- function(population,n.qtls=50,threshold=5,verbose=FALSE){
  if(missing(population)) stop("Please provide a population object\n")
  check.population(population)
  if(is.null(population$offspring$genotypes$real)){
    stop("No original genotypes in population$offspring$genotypes$real, load them in using intoPopulation function\n")
  }
  i <- 1
  qtls_found <- 0
  qtls <- NULL
  markers <- rownames(population$offspring$phenotypes)
  scores <- vector(mode="numeric",length=ncol(population$offspring$phenotypes))
  names(scores) <- colnames(population$offspring$phenotypes)
  while(qtls_found<n.qtls){
    cur_row <- matrix(QTLscan.internal(markers[i],t(population$offspring$phenotypes),population$offspring$genotypes$real),1,nrow(population$offspring$genotypes$real))
    cur_peaks <- getpeaks.internal(abs(cur_row),5)
    if(any(cur_peaks==2)){
      peakLocations <- which(cur_peaks==2)
      qtls_found <- qtls_found + length(peakLocations)
      if(verbose) cat(qtls_found,"qtls found\n")
      old_names<-names(qtls)
      qtls <- c(qtls,rep(i,length(peakLocations)))
      names(qtls) <- c(old_names,peakLocations)
    }
    i <- i+1
  }
  for(j in 1:qtls_found){
    group_a <- population$offspring$phenotypes[qtls[j],which(population$offspring$genotypes$real[as.numeric(names(qtls))[j],]==1)]
    group_b <- population$offspring$phenotypes[qtls[j],which(population$offspring$genotypes$real[as.numeric(names(qtls))[j],]==2)]
    scores <- scoreMixups.internal(group_a,scores,qtls_found)
    scores <- scoreMixups.internal(group_b,scores,qtls_found)
  }
  if(any(scores>threshold)){
    flagged <- which(scores>threshold)
    cat("Found",length(flagged),"possible mix-ups")
    for(flag in flagged){      
      cat(names(scores)[flag],":",scores[flag],"\n")
    }
  }
  invisible(scores)
}

scoreMixups.internal <- function(values,scores,qtls_found){
  cur_mean <- mean(values)
  cur_sd <- abs(sd(values))
  increase <- 1/qtls_found*100
  if(any(values>cur_mean+3*cur_sd)){
    positions <- names(values)[which(values>cur_mean+3*cur_sd)]
    scores[positions] <- scores[positions]+increase
  }
  if(any(values<cur_mean-3*cur_sd)){
    positions <- names(values)[which(values<cur_mean-3*cur_sd)]
    scores[positions] <- scores[positions]+increase
  }
  invisible(scores)
}
