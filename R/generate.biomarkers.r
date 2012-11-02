#
# findBiomarkers.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Mar, 2011
# Contains: findBiomarkers, pull.biomarkers, selectTopMarker.internal
#           scoreMarker.internal, convertfindBiomarkers.internal
#           splitPheno.internal, selectMarkersUsingMap.internal, 
#           filterGenotypes.internal, filterRow.internal, splitRowSubEM.internal
#

# findBiomarkers
#
# DESCRIPTION:
#  Function that chooses from the matrix only appropriate markers with specified rules
# PARAMETERS:
#   - population - Ril type object, must contain founders phenotypic data.
#   - orderUsing- which map should be used to order markers (default - none)
#     - map_genetic - genetic map
#     - map_physical - physical map
#   - treshold - If Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.
#   - overlapInd - Number of individuals that are allowed in the overlap
#   - proportion - Proportion of individuals expected to carrying a certain genotype 
#   - margin - Proportion is allowed to varry between this margin (2 sided)
#   - minChrLength - if maximal distance between the markers in the chromosome is lower than this value, whole chromosome will be dropped
#   - verbose - Be verbose
#   - debugMode - 1: Print our checks, 2: print additional time information
# OUTPUT:
#  An object of class cross
#
generate.biomarkers <- function(population, threshold=0.05, overlapInd = 10, proportion = c(50,50), margin = 15, verbose=FALSE, debugMode=0){
  if(missing(population)) stop("Population object not found.\n")
  check.population(population) # CHECK
  s<-proc.time()
  if(any(proportion < 1) || sum(proportion) != 100) stop("Wrong proportion paramete\n")
  if(overlapInd < 0 || overlapInd > ncol(population$offspring$phenotypes)) stop("overlapInd is a number (0,lenght of the row).")
  if(verbose && debugMode==1) cat("generate.biomarkers starting withour errors in checkpoint.\n")
  
  #*******CONVERTING CHILDREN PHENOTYPIC DATA TO GENOTYPES*******
  s1 <- proc.time()
  population <- generate.biomarkers.internal(population, threshold, overlapInd, proportion, margin, verbose, debugMode)
  e1 <- proc.time()
  if(verbose && debugMode==2)cat("Converting phenotypes to genotypes done in:",(e1-s1)[3],"seconds.\n")
  
  #*******RETURNING CROSS OBJECT*******
  e<-proc.time()
  if(verbose) cat("generate.biomarkers done in",(e-s)[3],"seconds\n")
  invisible(population)
}

############################################################################################################
#                  *** pull.biomarkers ***
#
# DESCRIPTION:
#  function returning all biomarkers or top marker matching given pattern
# 
# PARAMETERS:
#   population - an object of class population
#   pattern - vector of 0s and 1s (or 0,1,2s)
#   verbose - be verbose
# 
# OUTPUT:
#  vector/matrix
#
############################################################################################################
pull.biomarkers <- function(population,pattern,verbose=FALSE){
  if(missing(population)) stop("No population object provided.\n")
  if(is.null(population$offspring$genotypes$simulated)) stop("Population object doesn't contain de novo genotypes, run findBiomarkers.\n")
  markers <- population$offspring$genotypes$simulated
  if(verbose) cat("Selected",nrow(markers),"markers.\n")
  if(!missing(pattern)){
    if(length(pattern)!=ncol(markers)) stop("Wrong length of the pattern: ",length(pattern)," instead of: ",ncol(markers)," \n")
    if(verbose) cat("Selecting marker best matching given pattern.\n")
    markers <- selectTopMarker.internal(markers,pattern,verbose)
  }
  invisible(markers)
}

############################################################################################################
#                  *** selectTopMarker.internal  ***
#
# DESCRIPTION:
#  function returning all biomarkers or top marker matching given pattern
# 
# PARAMETERS:
#   population - an object of class population
#   pattern - vector of 0s and 1s (or 0,1,2s)
#   verbose - be verbose
# 
# OUTPUT:
#  vector/matrix
#
############################################################################################################
selectTopMarker.internal <- function(markers,pattern,verbose){
  markerPoints <- apply(markers,1,function(x){sum(x==pattern)})
  topMarker <- rownames(markers)[which.max(markerPoints)]
  if(verbose) cat("Markers best matching pattern:",topMarker,"with identity:",max(markerPoints)/ncol(markers)*100,"%\n")
  invisible(markers[topMarker,])
}

############################################################################################################
#                  *** convertfindBiomarkers.internal ***
#
# DESCRIPTION:
#  function splitting differentially expressed markers into two genotypes
# 
# PARAMETERS:
#   population - object of class population, must contain founders phenotypic data.
#   orderUsing- which map should be used to order markers (default - none)
#     - map_genetic - genetic map
#    - map_physical - physical map
#   treshold - if Rank Product pval for gene is lower that this value, we assume it is being diff. expressed.
#   overlapInd - number of individuals that are allowed in the overlap
#   proportion - proportion of individuals expected to carrying a certain genotype 
#   margin - proportion is allowed to varry between this margin (2 sided)
#   verbose - be verbose
#   debugMode - 1: Print our checks, 2: print additional time information 
# 
# OUTPUT:
#  object of class population
#
############################################################################################################
generate.biomarkers.internal <- function(population, treshold, overlapInd, proportion, margin, verbose=FALSE, debugMode=0){
  ### initialization
  populationType <- class(population)[2]
  if(verbose && debugMode==1) cat("generate.biomarkers.internal starting.\n")
  output <- NULL
  markerNames <- NULL 
  ### selection step
  ### up-regulated
  upNotNull <- which(population$founders$RP$pval[,1] > 0)
  upBelowTreshold <- which(population$founders$RP$pval[,1] < treshold)
  upSelected <- upBelowTreshold[which(upBelowTreshold%in%upNotNull)]
  print(length(upSelected))
  upParental <- population$founders$phenotypes[upSelected,]
  rownamesUp <- rownames(upParental)
  if(any(rownamesUp == "")) rownamesUp <- rownamesUp[-which(rownamesUp == "")]
  upRils <- population$offspring$phenotypes[rownamesUp,]
  ### down-regulated
  downNotNull <- which(population$founders$RP$pval[,2] > 0)
  downBelowTreshold <- which(population$founders$RP$pval[,2] < treshold)
  downSelected <- downBelowTreshold[which(downBelowTreshold%in%downNotNull)]
  print(length(downSelected))
  downParental <- population$founders$phenotypes[downSelected,]
  rownamesDown <- rownames(downParental)
  if(any(rownamesDown == "")) rownamesDown <- rownamesDown[-which(rownamesDown == "")]
  downRils <- population$offspring$phenotypes[rownamesDown,]
  
  ### checking if anything is selected and if yes - processing
  if(!(is.null(dim(upRils)))&&(nrow(upRils)!=0)){
    if(!(is.null(dim(downRils)))&&(nrow(downRils)!=0)){
      # best situation
      if(verbose) cat("Selected",nrow(downRils),"markers (DOWN), ",nrow(upRils),"markers (UP).\n")
      inupndown <- which(rownames(upRils) %in% rownames(downRils))
      if(verbose&&length(inupndown)>0){
        cat("WARNING: Overlap between UP n DOWN:",length(inupndown),", removing from UP.\n")
        upRils <- upRils[-inupndown,]
      }
      cur <- splitPheno.internal(downRils, downParental, overlapInd, proportion, margin, population$founders$groups, populationType, 0, 0, nrow(upRils),verbose)
      output <- rbind(output,cur)
    }else{
      if(verbose) cat("Selected ",nrow(upRils),"upregulated markers.\n")
    }
    cur <- splitPheno.internal(upRils, upParental, overlapInd, proportion, margin, population$founders$groups, populationType, 1, nrow(downRils), 0, verbose)
    output <- rbind(output,cur)
  }else{
    if(!(is.null(dim(downRils)))&&(nrow(downRils)!=0)){
      if(verbose) cat("Selected ",nrow(downRils),"downregulated markers.\n")
      cur <- splitPheno.internal(downRils, downParental, overlapInd, proportion, margin, population$founders$groups, populationType, 0, 0, 0,verbose)
      output <- rbind(output,cur)
    }else{
      stop("None of the markers was selected using specified treshold: ",treshold,"\n")
    }
  }
  
  ### putting results inside population object
  if(is.null(dim(output))) stop("No markers selected.")
  population$offspring$genotypes$simulated <- output
  colnames(population$offspring$genotypes$simulated) <- colnames(upRils)
  invisible(population)
}

############################################################################################################
#                  *** splitPheno.internal ***
#
# DESCRIPTION:
#  subfunction of convertfindBiomarkers.internal, splitting children markers using founders mean values
# 
# PARAMETERS:
#   offspring - matrix of up/down regulated genes in offspring
#   founders - matrix of up/down regulated genes in parents
#   overlapInd - Number of individuals that are allowed in the overlap
#   proportion - Proportion of individuals expected to carrying a certain genotype 
#   margin - Proportion is allowed to varry between this margin (2 sided)
#   groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
#   up - 1 - genes up 0 - down regulated
# 
# OUTPUT:
#  list containg genotype matrix and names of selected markers
#
############################################################################################################
#DANNY: TODO MERGE splitPhenoRowEM.internal into this function 
##K: left, I\'ll try to use apply here instead of for
splitPheno.internal <- function(offspring, founders, overlapInd, proportion, margin, groupLabels, populationType, up, done=0, left=0, verbose=FALSE){
  output <- NULL
  markerNames <- NULL
  s <-proc.time()
  for(x in 1:nrow(offspring)){
    cur <- splitPhenoRowEM.internal(x, offspring, founders, overlapInd, proportion, margin, groupLabels, up, populationType, verbose)
    if(!(is.null(cur))){
      output <- rbind(output,cur)
      markerNames <- c(markerNames,rownames(offspring)[x])
    }
    if(verbose){
      if((done+x)%%100==0){
        e <- proc.time()
        te <- ((e-s)[3]/x)*(nrow(offspring)-x+left)
        cat("Done with marker",done+x,"/",nrow(offspring)+left+done,". Time remaining:",te,"s\n")
      }
    }
  }
  rownames(output) <- markerNames
  invisible(output)
}

############################################################################################################
#                  *** splitPhenoRowEM.internal ***
#
# DESCRIPTION:
#  subfunction of splitRow.internal, splitting one row using EM algorithm
# 
# PARAMETERS:
#   x - name of currently processed row
#   offspring - matrix of up/down regulated genes in offspring
#   founders - matrix of up/down regulated genes in parents
#   overlapInd - Number of individuals that are allowed in the overlap
#   proportion - Proportion of individuals expected to carrying a certain genotype 
#   margin - Proportion is allowed to varry between this margin (2 sided)
#   groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
#   up - 1 - genes up 0 - down regulated
# 
# OUTPUT:
#  genotype row
#
############################################################################################################
splitPhenoRowEM.internal <- function(x, offspring, founders, overlapInd, proportion, margin, groupLabels, up=1, populationType, verbose=FALSE){
  y<-x
  x<-as.character(rownames(offspring)[x])
  aa <- tempfile()
  sink(aa)
  nrDistributions <- length(proportion)
  result <- rep(0,length(offspring[x,]))
  
  EM <- NULL
  s1<-proc.time()
  tryCatch(EM <- normalmixEM((offspring[x,]), k=nrDistributions, maxrestarts=0, maxit = 100,fast=FALSE),error = function(x){cat(x[[1]],"\n")})
  e1<-proc.time()
  sink()
  file.remove(aa)
  if(is.null(EM)){
        result <- NULL
  }else if(filterRow.internal(EM$lambda,proportion,margin)){
    if(populationType == "f2"){
      if(up==1){
        genotypes <- c(1:5)
      }else if(up==0){
        genotypes <- c(3,2,1,5,4)
      }
      for(i in (1:length(offspring[1,]))){
        if(any(EM$posterior[i,]>0.8)){
          result[i] <- genotypes[which.max(EM$posterior[i,])]
        }else if((EM$posterior[i,1]+EM$posterior[i,2])>0.8){
          result[i] <- genotypes[4]
        }else if((EM$posterior[i,2]+EM$posterior[i,3])>0.8){
          result[i] <- genotypes[5]
        }else{
          result[i] <- NA
        }
      }
    }else{
      if(up==1){
        genotypes <- c(1,2)
      }else if(up==0){
        genotypes <- c(2,1)
      }
      for(i in (1:length(offspring[1,]))){
        if(any(EM$posterior[i,]>0.8)){
          result[i] <- genotypes[which.max(EM$posterior[i,])]
        }else{
          result[i] <- NA
        }
      }
    }
    if(sum(is.na(result))>overlapInd){
      result <- NULL
    }
  }else{
   result <- NULL
  }
  invisible(result)
}


############################################################################################################
#                  *** filterRowSub.internal ***
#
# DESCRIPTION:
#   subfunction of filterGenotypes.internal, filtering one row
# 
# PARAMETERS:
#   genotypeRow - currently processed row
#   overlapInd - Number of individuals that are allowed in the overlap
#   proportion - Proportion of individuals expected to carrying a certain genotype 
#   margin - Proportion is allowed to varry between this margin (2 sided)
# 
# OUTPUT:
#  boolean
#
############################################################################################################
filterRow.internal<- function(lambda, proportion, margin){
  if(length(lambda)!=length(proportion)) return(FALSE)
  for(i in 1:length(lambda)){
    if((lambda[i]>((proportion[i]+margin/2)/100))||(lambda[i]<((proportion[i]-margin/2)/100))){
      return(FALSE)
    }
  }
  return(TRUE)
}
