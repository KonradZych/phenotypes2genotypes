############################################################################################################
#
# markerPlacementPlot.R
#
# Copyright (c) 2011, Danny Arends
#
# Modified by Konrad Zych
# 
# first written November 2011
# last modified December 2011
# last modified in version: 0.9.1
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
# Contains: markerPlacementPlot
#
############################################################################################################

markerPlacementPlot <- function(population, placeUsing=c("qtl","correlation"), map=c("genetic","physical"), cross){
  if(missing(population)) stop("Please provide a population object\n")
  if(is.null(population$offspring$genotypes$real)){
    stop("No original genotypes in population$offspring$genotypes$real, load them in using intoPopulation\n")
  }
  check.population(population)
  if(missing(cross)){
    cross <- genotypesToCross.internal(population, "simulated")
  }
  if(placeUsing=="correlation"){
    cormatrix <- map2mapCorrelationMatrix(cross,population)
    s <- NULL
    p <- seq(1,5,0.1)
    for(x in p){
      maximums <- apply(abs(cormatrix),2,max)
      means <- apply(abs(cormatrix),2,mean)
      sds <- apply(abs(cormatrix),2,sd)
      selected <- which(maximums > (means+x*sds))
      s <- c(s,length(selected))
    }
  plot(p,s,type='o',main="Number of markers placed",xlab="corThreshold",ylab="# of markers")
  }else{
    map <- checkParameters.internal(map,c("genetic","physical"),"map")
    if(map=="genetic"){
      cur_map <- population$maps$genetic
    }else{
      cur_map <- population$maps$physical
    }
    genotypes <- population$offspring$genotypes$real
    markers <- markernames(cross)
    phenotypes <- pull.geno(cross)[,markers]
    results <- NULL
    rnames <- NULL
    for(x in 1:ncol(phenotypes)){
      if(x%%10==0) cat("-----",x,"-------\n")
      results <- rbind(results,apply(genotypes,1, 
        function(geno){
          linmod <- lm(phenotypes[,x] ~ geno)
          -log10(anova(linmod)[[5]][1])
        }
      ))
      rnames <- c(rnames,colnames(phenotypes)[x])
    }
  rownames(results) <- rnames
  singleqtl <- NULL
  noqtl <- NULL
  multipleqtl <- NULL
  thrRange <- seq(0,20,0.1)
  for(tr in thrRange){
    peaks <- getpeaks.internal(results,tr)
    rownames(peaks) <- rownames(results)
    colnames(peaks) <- colnames(results)
    nqtls <- checkpeaks.internal(results,cur_map,tr)
    singleqtl <- c(singleqtl,sum(nqtls==1))
    noqtl <- c(noqtl,sum(nqtls==0))
    multipleqtl <- c(multipleqtl,sum(nqtls>1))
  }
  plot(thrRange,noqtl,type='o',col="red",main="Number of markers placed",xlab="threshold",ylab="# of markers",ylim=c(0,ncol(phenotypes)))
  points(thrRange,singleqtl,type='o',col="green")
  points(thrRange,multipleqtl,type='o',col="blue")
  legend(c(10,18),c(0.85*ncol(phenotypes),ncol(phenotypes)),c("no peak","single peak","multiple peaks"),col=c("red","green","blue"),cex=0.8,pch=21,lwd=2,bg="white")
  invisible(results)
  }
}

getpeaks.internal <- function(qtlprofiles, cutoff = 4.0){
  cat("Starting peak detection above",cutoff,"\n")
  mmatrix <- NULL
  #qtlprofiles[which(qtlprofiles==Inf)]<-10000
  for(x in 1:nrow(qtlprofiles)){
    peak <- FALSE
    curmax <- 0
    curmaxindex <- 1
    marker <- 1
    maximums <- NULL
    mrow <- rep(0,ncol(qtlprofiles))
    for(ab in (qtlprofiles[x,]>cutoff | qtlprofiles[x,]<(-cutoff))){
      if(ab){
        peak <- TRUE
        if(qtlprofiles[x,marker]/abs(qtlprofiles[x,marker]) > 0){
          if(qtlprofiles[x,marker] > curmax){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- marker
          }
        }else{
          if(qtlprofiles[x,marker] < (-curmax)){
            curmax <- qtlprofiles[x,marker]
            curmaxindex <- -marker
          }
        }
        if(ncol(qtlprofiles)==marker){
          if(curmax!=0) maximums <- c(maximums,curmaxindex)
        }
      }else{
        if(curmax!=0) maximums <- c(maximums,curmaxindex)
        peak <- FALSE
        curmax <- 0
      }
      marker <- marker+1
    }
    mrow[which(qtlprofiles[x,] > cutoff)] <- 1
    mrow[which(qtlprofiles[x,] < -cutoff)] <- -1
    for(a in which(maximums>0)){
      mrow[maximums[a]] <- 2
    }
    for(b in which(maximums<0)){
      mrow[(-maximums[b])] <- -2
    }
    mmatrix <- rbind(mmatrix,mrow)
  }
  mmatrix
}

checkpeaks.internal <- function(results,cur_map,threshold){
  result <- NULL
  for(i in 1:nrow(results)){
  markerpeaks <- NULL
    for(chr in unique(cur_map[,1])){
      markers <- rownames(cur_map)[which(cur_map[,1]==chr)]
      markerpeaks <- c(markerpeaks,(sum(results[i,markers]>threshold)/length(markers)))
    }
  result <- rbind(result,markerpeaks)
  }
  invisible(result)
}