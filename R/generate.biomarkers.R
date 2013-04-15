#
# generate.biomarkers.R
#
# Copyright (c) 2010-2013 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified April, 2013
# first written Mar, 2011
# Contains: generate.biomarkers, pull.biomarkers, selectTopMarker.internal
#           scoreMarker.internal, convertfindBiomarkers.internal
#           splitPheno.internal, selectMarkersUsingMap.internal, 
#           filterGenotypes.internal, filterRow.internal, splitRowSubEM.internal
#

# generate.biomarkers
#
# DESCRIPTION:
#  Function that chooses from the matrix only appropriate markers with specified rules
# PARAMETERS:
#   - population - An object of class population.
#   - threshold - If  pval for gene is lower that this value, we assume it is being diff. expressed.
#   - overlapInd - Number of individuals that are allowed in the overlap.
#   - proportion - Proportion of individuals expected to carrying a certain genotype. 
#   - pProb - Posterior probability threshold used to assugn genotypes.
#   - env - Vector contatining information about environment for each of the individuals in the set (numeric value).
#   - verbose - Be verbose.
#   - debugMode - 1: Print our checks, 2: print additional time information
# OUTPUT:
#  An object of class cross
#
generate.biomarkers <- function(population, threshold=0.05, overlapInd = 10, proportion = c(50,50), margin = 15, pProb=0.8, env, verbose=FALSE, debugMode=0){
  
  ### checks
  ### check population
  if(missing(population)) stop("Population object not found.\n")
  check.population(population)

  ### check threshold
  if(!is.numeric(threshold))                        stop("threshold should be a numeric value\n")
  if(threshold < 0)                                 stop("threshold should be a positivite value\n")

  ### check overlapInd
  if(!is.numeric(overlapInd))                       stop("overlapInd should be a numeric value\n")
  if(overlapInd < 0)                                stop("overlapInd should be a positivite value\n")

  ### check proportion
  if(any(!is.numeric(proportion)))                  stop("overlapInd should be a numeric value\n")
  if(any(proportion < 1) || sum(proportion) != 100) stop("Proportion should be > 0 and < 100 and sum up to 100\n")

  ### check margin
  if(!is.numeric(margin))                           stop("margin should be a numeric value\n")
  if(margin < 0 || margin > 100)                    stop("margin should be a number between 0 and 100\n")
  
  ### check pProb
  if(!is.numeric(pProb))                            stop("pProb should be a numeric value\n")
  if(pProb < 0 || pProb > 1)                        stop("pProb should be a number between 0 and 1\n")

  if(verbose && debugMode==1) cat("generate.biomarkers starting withour errors in checkpoint.\n")
  s  <- proc.time()
  
  #*******CONVERTING CHILDREN PHENOTYPIC DATA TO GENOTYPES*******
  s1 <- proc.time()
  if(!is.null(population$annots)){ #TODO: Merge this 2 functions into 1
    population <- generate.biomarkers.onthefly.internal(population, threshold, overlapInd, proportion, margin, pProb, verbose, debugMode)
  }else{
    population <- generate.biomarkers.internal(population, threshold, overlapInd, proportion, margin, pProb, verbose, debugMode)
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
generate.biomarkers.internal <- function(population, treshold, overlapInd, proportion, margin, pProb=0.8, env, verbose=FALSE, debugMode=0){
  ### initialization
  populationType <- class(population)[2]
  if(verbose && debugMode==1) cat("generate.biomarkers.internal starting.\n")
  output      <- NULL
  outputEM    <- vector(2,mode="list")
  markerNames <- NULL 
  
  ### selection step
  ### checking if any of the phenotypes is up-regulated
  upRegulatedPhenos      <- selectPhenotypes(population, treshold, RPcolumn=1)
  
  ### checking if any of the phenotypes is down-regulated
  downRegulatedPhenos    <- selectPhenotypes(population, treshold, RPcolumn=2)
  
  if(!is.null(rownames(upRegulatedPhenos)) && !is.null(rownames(downRegulatedPhenos))){
    inupndown                    <- which(rownames(upRegulatedPhenos) %in% rownames(downRegulatedPhenos))
    if(length(inupndown)>0)      upRegulatedPhenos      <- upRegulatedPhenos[-inupndown,]
  }

  ### if any of the phenotypes is up-regulated - process them
  if(!(is.null(dim(upRegulatedPhenos)))&&(nrow(upRegulatedPhenos)!=0)){
    if(verbose) cat("Selected ",nrow(upRegulatedPhenos),"upregulated potential markers.\n")
    cur                     <- splitPheno.internal(upRegulatedPhenos, upParental, overlapInd, proportion, margin, pProb, populationType, 1, 0, 0, verbose)
    output                  <- rbind(output,cur[[1]])
    outputEM[[2]]           <- cur[[2]]
    names(outputEM[[2]])    <- rownames(population$offspring$phenotypes)
    
  }else{
    if(verbose) cat("Selected none upregulated potential markers.\n")
  }

  ### if any of the phenotypes is down-regulated - process them
  if(!(is.null(dim(downRegulatedPhenos)))&&(nrow(downRegulatedPhenos)!=0)){
    if(verbose) cat("Selected ",nrow(downRegulatedPhenos),"downregulated potential markers.\n")
    cur                     <- splitPheno.internal(downRegulatedPhenos, downParental, overlapInd, proportion, margin, pProb, populationType, 0, 0, 0,verbose)
    output                  <- rbind(output,cur[[1]])
    outputEM[[1]]           <- cur[[2]]
    names(outputEM[[1]])    <- rownames(population$offspring$phenotypes)
  }else{
    if(verbose) cat("Selected none downregulated potential markers.\n")
  }
  
  if(verbose) cat("Generated ",nrow(output),"markers.\n")
  ### putting results inside population object
  if(is.null(dim(output))) stop("No markers selected.")

  population$offspring$genotypes$simulated <- output
  population$offspring$genotypes$EM <- outputEM
  colnames(population$offspring$genotypes$simulated) <- colnames(upRegulatedPhenos)
  invisible(population)
}

# selectPhenotypes
#
# DESCRIPTION:
#  Function that selects offsprings fulfilling the criteria
# PARAMETERS:
#   - population - An object of class population.
#   - threshold - If  pval for gene is lower that this value, we assume it is being diff. expressed.

# OUTPUT:
#  A matrix with selected phenotypes
#
selectPhenotypes <- function(population, treshold, RPcolumn){
  notNullPhenotypes   <- which(population$founders$RP$pval[,RPcolumn] > 0)      # rank product gives a score for 0 sometimes -> this is below the threshold but these phenotypes wshould not be selected
  belowTreshold       <- which(population$founders$RP$pval[,RPcolumn] < treshold) # phenos diff expressed with pval lower than threshold
  selected            <- belowTreshold[which(belowTreshold%in%notNullPhenotypes)]
  selectedParental    <- population$founders$phenotypes[selected,]
  rownamesOfSelected  <- rownames(selectedParental)
  if(any(rownamesOfSelected == "")) rownamesOfSelected <- rownamesOfSelected[-which(rownamesOfSelected == "")]
  selectedRils              <- population$offspring$phenotypes[rownamesOfSelected,]
  invisible(selectedRils)
}

t.test_line.internal <- function(dataRow,threshold){
  idx <- which(!is.na(dataRow))
  y <- dataRow[idx]
  half <- floor(length(y)/2)
  res <- t.test(sort(y)[1:half],sort(y)[half:length(y)])
  if(res$p.value<threshold) invisible(TRUE)
  invisible(FALSE)
}

#TODO: This function is equal to check.and.generate.internal merge them !!!
generate.internal <- function(dataRow, pval, threshold, overlapInd, proportion, margin, pProb, curlineNR, populationType, verbose){
  cur <- NULL
  dataRow <- as.numeric(dataRow)
  if(pval[1]<threshold && pval[1]!=0){
    cur <- splitPhenoRowEM.internal(dataRow, overlapInd, proportion, margin, pProb, 0, populationType, verbose)[[1]]
  }else if(pval[2]<threshold && pval[2]!=0){
    cur <- splitPhenoRowEM.internal(dataRow, overlapInd, proportion, margin, pProb, 1, populationType, verbose)[[1]]
  }
  if(!is.null(cur)) cur <- c(curlineNR,cur)

  invisible(cur)
}

#TODO: This function is equal to generate.internal merge them !!!
check.and.generate.internal <- function(dataRow, threshold, overlapInd, proportion, margin, pProb, curlineNR, populationType, verbose){
  cur <- NULL
  dataRow <- as.numeric(dataRow)
  if(t.test_line.internal(dataRow,threshold)){
     cur <- splitPhenoRowEM.internal(dataRow, overlapInd, proportion, margin, pProb, 1, populationType, verbose)[[1]]
     if(!is.null(cur)) cur <- c(curlineNR,cur)
  }
  invisible(cur)
}

generate.biomarkers.onthefly.internal <- function(population, threshold, overlapInd, proportion, margin, pProb=0.8, env, verbose=FALSE, debugMode=0){
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
          genoLine <- generate.internal(analysedLines[[analysedLineNR]][-1], population$founders$RP$pval[probeID,], threshold, overlapInd, proportion, margin, pProb, curlineNR, populationType, verbose)
        }else{
         genoLine <- check.and.generate.internal(analysedLines[[analysedLineNR]][-1], threshold, overlapInd, proportion, margin, pProb, curlineNR, populationType, verbose)
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
            }else{
              #TODO: What do we need to do when we are in the ELSE ?
            }
            newProbe[idb] <- probeB_[idb]
            probe_ <- newProbe
          }
        }else{
          #TODO: What do we need to do when we are in the ELSE ?
        }
      }else{
        #TODO: What do we need to do when we are in the ELSE ?
      }
      newGeno <- rbind(newGeno,c(probeName,probe_))
    }
  }else{
    #TODO: What do we need to do when we are in the ELSE ?
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
splitPheno.internal <- function(offspring, founders, overlapInd, proportion, margin, pProb=0.8, populationType, up, done=0, left=0, verbose=FALSE){
  output        <- NULL
  outputEM      <- NULL
  markerNames   <- NULL
  s             <- proc.time()
  printedProc   <- NULL
  for(x in 1:nrow(offspring)){
    cur <- splitPhenoRowEM.internal(offspring[x,], overlapInd, proportion, margin, pProb, up, populationType, verbose)
    if(!(is.null(cur[[1]]))){
      output      <- rbind(output,cur[[1]])
      markerNames <- c(markerNames,rownames(offspring)[x])
    }
    outputEM    <- rbind(outputEM,cur[[2]])
    if(verbose){
      perc <- round((x+done)*100/nrow(offspring+done+left))
      if(perc%%10==0 && !(perc%in%printedProc)){
        e <- proc.time()
        printedProc <- c(printedProc,perc)
        te <- ((e-s)[3]/x)*(nrow(offspring)-x+left)
        cat("Done with: ",perc,"%. Estimated time remaining:",te,"s\n")
      }
    }
  }
  rownames(output) <- markerNames

  invisible(list(output,outputEM))
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
splitPhenoRowEM.internal <- function(x, overlapInd, proportion, margin, pProb=0.8, up=1, populationType, verbose=FALSE){
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
        if(any(EM$posterior[i,]>pProb)){
          result[idx[i]] <- genotypes[which.max(EM$posterior[i,])]
        }else if((EM$posterior[i,1]+EM$posterior[i,2])>pProb){
          result[idx[i]] <- genotypes[4]
        }else if((EM$posterior[i,2]+EM$posterior[i,3])>pProb){
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
        if(any(EM$posterior[i,]>pProb)) result[idx[i]] <- genotypes[which.max(EM$posterior[i,])]
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

