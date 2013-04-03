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
generate.biomarkers <- function(population, threshold=0.05, overlapInd = 10, proportion = c(50,50), margin = 15, p.prob=0.8, env, verbose=FALSE, debugMode=0){
  if(missing(population)) stop("Population object not found.\n")
  check.population(population) # CHECK
  s<-proc.time()
  if(any(proportion < 1) || sum(proportion) != 100) stop("Proportion should be > 0 and < 100")

  if(verbose && debugMode==1) cat("generate.biomarkers starting withour errors in checkpoint.\n")
  
  #*******CONVERTING CHILDREN PHENOTYPIC DATA TO GENOTYPES*******
  s1 <- proc.time()
  if(!is.null(population$annots)){ #TODO: Merge this 2 functions into 1
    population <- generate.biomarkers.onthefly.internal(population, threshold, overlapInd, proportion, margin, p.prob, verbose, debugMode)
  }else{
    population <- generate.biomarkers.internal(population, threshold, overlapInd, proportion, margin, p.prob, verbose, debugMode)
  }
  e1 <- proc.time()
  if(verbose && debugMode==2) cat("Converting phenotypes to genotypes done in:",(e1-s1)[3],"seconds.\n")
  
  #*******RETURNING CROSS OBJECT*******
  e<-proc.time()
  if(verbose) cat("generate.biomarkers done in",(e-s)[3],"seconds.\n")
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
generate.biomarkers.internal <- function(population, treshold, overlapInd, proportion, margin, p.prob=0.8, env, verbose=FALSE, debugMode=0){
  ### initialization
  populationType <- class(population)[2]
  if(verbose && debugMode==1) cat("generate.biomarkers.internal starting.\n")
  output <- NULL
  outputEM <- vector(2,mode="list")
  markerNames <- NULL 
  ### selection step
  ### up-regulated
  upNotNull <- which(population$founders$RP$pval[,1] > 0)
  upBelowTreshold <- which(population$founders$RP$pval[,1] < treshold)
  upSelected <- upBelowTreshold[which(upBelowTreshold%in%upNotNull)]
  upParental <- population$founders$phenotypes[upSelected,]
  rownamesUp <- rownames(upParental)
  if(any(rownamesUp == "")) rownamesUp <- rownamesUp[-which(rownamesUp == "")]
  upRils <- population$offspring$phenotypes[rownamesUp,]
  ### down-regulated
  downNotNull <- which(population$founders$RP$pval[,2] > 0)
  downBelowTreshold <- which(population$founders$RP$pval[,2] < treshold)
  downSelected <- downBelowTreshold[which(downBelowTreshold%in%downNotNull)]
  downParental <- population$founders$phenotypes[downSelected,]
  rownamesDown <- rownames(downParental)
  if(any(rownamesDown == "")) rownamesDown <- rownamesDown[-which(rownamesDown == "")]
  downRils <- population$offspring$phenotypes[rownamesDown,]
  
  ### checking if anything is selected and if yes - processing
  if(!(is.null(dim(upRils)))&&(nrow(upRils)!=0)){
    if(!(is.null(dim(downRils)))&&(nrow(downRils)!=0)){
      # best situation
      if(verbose) cat("Selected",nrow(downRils),"potential markers (downregulated), ",nrow(upRils),"potential markers (upregulated).\n")
      inupndown <- which(rownames(upRils) %in% rownames(downRils))
      if(verbose&&length(inupndown)>0){
        #cat("WARNING: Overlap between UP n DOWN:",length(inupndown),", removing from UP.\n")
        upRils <- upRils[-inupndown,]
      }
      cur <- splitPheno.internal(downRils, downParental, overlapInd, proportion, margin, p.prob, populationType, 0, 0, nrow(upRils),verbose)
      output <- rbind(output,cur)#[[1]])
      #outputEM[[1]] <- cur[[2]]
     # names(outputEM[[1]]) <- rownames(downRils)
    }else{
      if(verbose) cat("Selected ",nrow(upRils),"upregulated potential markers.\n")
    }
    cur <- splitPheno.internal(upRils, upParental, overlapInd, proportion, margin, p.prob, populationType, 1, 0, 0, verbose)
    output <- rbind(output,cur)#[[1]])
    #outputEM[[2]] <- cur[[2]]
    #names(outputEM[[2]]) <- rownames(upRils)
  }else{
    if(!(is.null(dim(downRils)))&&(nrow(downRils)!=0)){
      if(verbose) cat("Selected ",nrow(downRils),"downregulated potential markers.\n")
      cur <- splitPheno.internal(downRils, downParental, overlapInd, proportion, margin, p.prob, populationType, 0, 0, 0,verbose)
      output <- rbind(output,cur)#[[1]])
      #outputEM[[1]] <- cur[[2]]
      #names(outputEM[[1]]) <- rownames(downRils)
    }else{
      stop("None of the markers was selected using specified treshold: ",treshold,"\n")
    }
  }
  if(verbose) cat("Generated ",nrow(output),"markers.\n")
  ### putting results inside population object
  if(is.null(dim(output))) stop("No markers selected.")
 # if(length(table(env))>1){
  #  print("multiple")
   # if(populationType=="f2"){
   #   output <- t(apply(output,1,cleanGeno.internal,env,c(1,2,3)))
   # }else{
   #   output <- t(apply(output,1,cleanGeno.internal,env,c(1,2)))
   # }
 # }
  population$offspring$genotypes$simulated <- output
  #population$offspring$genotypes$EM <- outputEM
  colnames(population$offspring$genotypes$simulated) <- colnames(upRils)
  invisible(population)
}

t.test_line.internal <- function(dataRow,threshold){
  idx <- which(!is.na(dataRow))
  y <- dataRow[idx]
  half <- floor(length(y)/2)
  res <- t.test(sort(y)[1:half],sort(y)[half:length(y)])
  if(res$p.value<threshold){
    invisible(TRUE)
  }else{
    invisible(FALSE)
  }
}

generate.internal <- function(dataRow, pval, threshold, overlapInd, proportion, margin, p.prob, curlineNR, populationType, verbose){
  cur <- NULL
  dataRow <- as.numeric(dataRow)
  if(pval[1]<threshold && pval[1]!=0){
    cur <- splitPhenoRowEM.internal(dataRow, overlapInd, proportion, margin, p.prob, 0, populationType, verbose)[[1]]
  }else if(pval[2]<threshold && pval[2]!=0){
    cur <- splitPhenoRowEM.internal(dataRow, overlapInd, proportion, margin, p.prob, 1, populationType, verbose)[[1]]
  }
  if(!is.null(cur)){
    cur <- c(curlineNR,cur)
  }
  invisible(cur)
}

check.and.generate.internal <- function(dataRow, threshold, overlapInd, proportion, margin, p.prob, curlineNR, populationType, verbose){
  cur <- NULL
  dataRow <- as.numeric(dataRow)
  if(t.test_line.internal(dataRow,threshold)){
     cur <- splitPhenoRowEM.internal(dataRow, overlapInd, proportion, margin, p.prob, 1, populationType, verbose)[[1]]
     if(!is.null(cur)){
       cur <- c(curlineNR,cur)
     }
  }
  invisible(cur)
}

generate.biomarkers.onthefly.internal <- function(population, threshold, overlapInd, proportion, margin, p.prob=0.8, env, verbose=FALSE, debugMode=0){
  analysedFile <- file(population$offspring$phenotypes,"r")
  header <- unlist(strsplit(readLines(analysedFile,n=1),"\t"))
  genoMatrix <- NULL
  phenoMatrix <- NULL
  lineNR <- 0
  populationType <- class(population)[2]
  analysedLines <- readLines(analysedFile,n=population$sliceSize)
  while(!((length(analysedLines) == 0) && (typeof(analysedLines) == "character"))){ #TODO: Wrong condition in the while
    analysedLines <- strsplit(analysedLines,"\t")
    eo <- proc.time()
    if(!((length(analysedLines) == 0) && (typeof(analysedLines) == "character"))){  #TODO: Wrong condition in the if (should not be here)
     cat("------",length(analysedLines),"\n")
      for(analysedLineNR in 1:length(analysedLines)){
        curlineNR <- analysedLineNR + lineNR
        if(curlineNR %% 500 == 0){
          en <- proc.time()
          cat("Analysing phenotype:",curlineNR,"time:",(en-eo)[3],"s\n")
        }
        probeID <- as.character(population$annots[curlineNR,1])
        if(!(is.null(population$founders$RP$pval[probeID,]))){
          genoLine <- generate.internal(analysedLines[[analysedLineNR]][-1], population$founders$RP$pval[probeID,], threshold, overlapInd, proportion, margin, p.prob, curlineNR, populationType, verbose)
        }else{
         genoLine <- check.and.generate.internal(analysedLines[[analysedLineNR]][-1], threshold, overlapInd, proportion, margin, p.prob, curlineNR, populationType, verbose)
        }
        if(!is.null(genoLine)){
          genoMatrix <- rbind(genoMatrix,genoLine)
          phenoMatrix <- rbind(phenoMatrix,c(analysedLineNR,analysedLines[[analysedLineNR]][-1]))
        }
      }
      
      lineNR <- lineNR+length(analysedLines)
      analysedLines <- readLines(analysedFile,n=population$sliceSize)
    }
  }
  genoMatrix <- mergeEnv.internal(population, genoMatrix)
  rownames(genoMatrix) <- genoMatrix[,1]
  rownames(phenoMatrix) <- phenoMatrix[,1]
  genoMatrix <- genoMatrix[,-1]
  phenoMatrix <- phenoMatrix[,-1]
  colnames(genoMatrix) <- header
  colnames(phenoMatrix) <- header
  population$offspring$phenotypes <- phenoMatrix
  population$offspring$genotypes$simulated <- genoMatrix
  close(analysedFile)
  invisible(population)
}

mergeEnv.internal <- function(population, genoMatrix){
  ### check if there is anything to merge and if so -> merge
  if(length(unique(population$annots[,3]))<nrow(population$annots)){
    done <- NULL
    newGeno <- NULL
    for(probenr in 1:nrow(genoMatrix)){
      probe <- genoMatrix[probenr,]
      probeID <- population$annots[probe[1],2]
      probeName <- population$annots[probe[1],1]
      probe_ <- probe[-1]
      if(!(probeID %in% done)){
        done <- c(done,probeID)
        probes <- which(population$annots[,2]==probeID)
        cat(probes,":",length(probes),"\n")
        if(length(probes)>1){
          newProbe <- probe
          idx <- which(!is.na(probe_))
          for(probeB in probes[2:length(probes)]){
            cat("Mergining:",probe[1],"with",probeB[1],"\n")
            probeB_ <- probeB[-1]
            idb <- which(!(is.na(probeB_)))
            if(any(idb%in%idx)){
              newProbe[which(idb%in%idx)] <- mean(newProbe[which(idb%in%idx)], probeB_[which(idb%in%idx)])
              idb <- idb[which(!(idb%in%idx))]
            }
            newProbe[idb] <- probeB_[idb]
            probe_ <- newProbe
          }
        }
      }
      newGeno <- rbind(newGeno,c(probeName,probe_))
    }
  }
  invisible(newGeno)
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
#DANNY: TODO MERGE splitPhenoRowEM.internal into this function  <- TODO: Do somehting with the TODOs
##K: left, I\'ll try to use apply here instead of for
splitPheno.internal <- function(offspring, founders, overlapInd, proportion, margin, p.prob=0.8, populationType, up, done=0, left=0, verbose=FALSE){
  output      <- NULL
  markerNames <- NULL
  s           <- proc.time()
  for(x in 1:nrow(offspring)){
    cur <- splitPhenoRowEM.internal(offspring[x,], overlapInd, proportion, margin, p.prob, up, populationType, verbose)
    if(!(is.null(cur[[1]]))){
      output <- rbind(output,cur[[1]])
      markerNames <- c(markerNames,rownames(offspring)[x])
    }

    if(verbose){
      perc <- round((x+done)*100/nrow(offspring+done+left))
      if(perc%%10==0){
        e <- proc.time()
        te <- ((e-s)[3]/x)*(nrow(offspring)-x+left)
        cat("Done with: ",perc,"%. Estimated time remaining:",te,"s\n")
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
splitPhenoRowEM.internal <- function(x, overlapInd, proportion, margin, p.prob=0.8, up=1, populationType, verbose=FALSE){
  #cat("post.prob = ",p.prob,"\n")
  aa <- tempfile()
  sink(aa)                                   #TODO: When we sink, we need to try{}catch so that we can undo our sink even when an error occurs
  nrDistributions <- length(proportion)
  result <- rep(NA,length(x))
  
  EM <- NULL
  idx <- which(!(is.na(x)))
  idw <- length(which((is.na(x))))
  y <- x[idx]
  s1<-proc.time()
  print(tryCatch)
  tryCatch(EM <- normalmixEM(y, k=nrDistributions, lambda= proportion, maxrestarts=1, maxit = 300, fast=FALSE),error = function(x){cat(x[[1]],"\n")})
  e1<-proc.time()
  sink()
  file.remove(aa)
  result <- NULL
  if(filterRow.internal(EM$lambda,proportion,margin)){
    if(populationType == "f2"){
      genotypes <- c(1:5)
      if(up==0) genotypes <- c(3,2,1,5,4)
      for(i in (1:length(y))){
        if(any(EM$posterior[i,]>p.prob)){
          result[idx[i]] <- genotypes[which.max(EM$posterior[i,])]
        }else if((EM$posterior[i,1]+EM$posterior[i,2])>p.prob){
          result[idx[i]] <- genotypes[4]
        }else if((EM$posterior[i,2]+EM$posterior[i,3])>p.prob){
          result[idx[i]] <- genotypes[5]
        }else{
          result[idx[i]] <- NA
        }
      }
    }else{
      genotypes <- c(1:2)
      if(up==0) genotypes <- c(2,1)
      for(i in (1:length(y))){
        result[idx[i]] <- NA
        if(any(EM$posterior[i,]>p.prob)) result[idx[i]] <- genotypes[which.max(EM$posterior[i,])]
      }
    }
    if((sum(is.na(result))-idw)>overlapInd) result <- NULL
  }
  invisible(list(result,EM))
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
  if(length(lambda) != length(proportion)) return(FALSE)
  for(i in 1:length(lambda)){
    if((lambda[i]>((proportion[i]+margin/2)/100)) || (lambda[i]<((proportion[i]-margin/2)/100))) return(FALSE)    #TODO use abs and merge the 2 similar IF conditions
  }
  return(TRUE)
}

