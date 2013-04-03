#
# cross.saturate.R
#
# Copyright (c) 2010-2013 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified January, 2013
# first written Mar, 2011
# Contains: cross.saturate, rearrangeMarkers, bestCorelated.internal
#           map2mapCorrelationMatrix.internal, map2mapImage
#

# cross.saturate
#
# DESCRIPTION:
#  Saturate an existing genetic map by adding markers derived from expression
# OUTPUT:
#  An object of class cross
#
cross.saturate <- function(population, cross, map=c("genetic","physical"), placeUsing=c("qtl","correlation"), flagged = c("remove","warn","ignore"), model, threshold=3, chr, env, keep.redundant=FALSE, use.orderMarkers=FALSE, verbose=FALSE, debugMode=0){
  if(missing(population)) stop("Please provide a population object\n")
  flagged <- match.arg(flagged)
  if(missing(env)) env <- rep(1,ncol(population$offspring$phenotypes))
  if(length(env)!=ncol(population$offspring$phenotypes)) stop("Incorrect environmental vector!\n")
  populationType <- class(population)[2]
  check.population(population)
  if(!is.numeric(threshold)||is.na(threshold)) stop("Please provide correct threshold")
  if(threshold<0) stop("Threshold needs to be > 0")
  if(placeUsing=="correlation" && threshold >= 5) cat("WARNING: threshold too high, few new markers will be selected\n")
  if(placeUsing=="qtl" && threshold >= 20) cat("WARNING: threshold too high, few new markers will be selected\n")

  map <- match.arg(map)
  placeUsing <- checkParameters.internal(placeUsing,c("qtl","correlation"),"placeUsing")
  if(missing(cross)){
    if(is.null(population$offspring$genotypes$real)) stop("No original genotypes in population$offspring$genotypes$real, load them in using add.to.population")
    if(is.null(population$offspring$genotypes$simulated)) stop("No genotype data in population$offspring$genotypes$simulated, run generate.biomarkers first")
    #*******SAVING CROSS OBJECT*******
    s1 <- proc.time()
    aa <- tempfile()
    sink(aa)
    cross <- genotypesToCross.internal(population,"simulated",verbose=verbose,debugMode=debugMode)
    sink()
    file.remove(aa)
    e1 <- proc.time()
    if(verbose && debugMode==2)cat("Saving data into cross object done in:",(e1-s1)[3],"seconds.\n")
  }else{
    population <- set.geno.from.cross(cross,population,map)
    population <- scan.qtls(population,map,env=env)
    aa <- tempfile()
    sink(aa)
    cross <- genotypesToCross.internal(population,"simulated",verbose=verbose,debugMode=debugMode)
    sink()
    file.remove(aa)
    e1 <- proc.time()
  }
  if(!(all(rownames(population$offspring$genotypes$simulated)%in%rownames(population$offspring$genotypes$qtl$lod)))) stop("QTL scan results don't match with simulated genotypes, please, run scan.qtls function")
  if(!(all(rownames(population$offspring$genotypes$qtl$lod)%in%rownames(population$offspring$genotypes$simulated)))) stop("QTL scan results don't match with simulated genotypes, please, run scan.qtls function")

 if(map=="genetic"){
    matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$genetic))
    if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
    if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
      population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
      population$maps$genetic <- population$maps$genetic[rownames(population$offspring$genotypes$real),]
      if(verbose) cat(nrow(population$offspring$genotypes$real)-length(matchingMarkers),"markers were removed due to name mismatch\n")
    }
    cur_map <- population$maps$genetic
  }else{
    matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$physical))
    if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
    if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
      population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
      population$maps$physical <- population$maps$physical[rownames(population$offspring$genotypes$real),]
      if(verbose) cat(nrow(population$offspring$genotypes$real)-length(matchingMarkers),"markers were removed due to name mismatch\n")
    }
    cur_map <- population$maps$physical
  }
  n.originalM <- nrow(population$offspring$genotypes$real)
  ### saturating only a subset of chromosomes
  if(missing(chr)){
    if(verbose) cat("Saturating all the chromosomes in the set\n")
    chr = unique(cur_map[,1])
  }else{
    availableChr = unique(cur_map[,1])
    if(any(!(chr%in%availableChr))) stop("Incorrect chr parameter!\n")
    if(verbose) cat("Saturating chromosomes:\n",paste(chr,",",sep=""),"\n")
  }
  #*******ENRICHING ORIGINAL MAP*******
  s1 <- proc.time()
  cross <- rearrangeMarkers(cross, population, populationType, cur_map, threshold, placeUsing, flagged, env, addMarkers=TRUE, keep.redundant, chr, verbose=verbose)
  count <- 1
  while(cross$left>1000){ # TODO: What does this 1000 do ??? Please explain and make it a user defined parameter
    count <- count+1      # TODO: We dont start counting from 2 ???!!! we start at 0
    population <- scan.qtls(set.geno.from.cross(cross, population, map),map)
    aa <- tempfile()
    sink(aa)
    cross <- genotypesToCross.internal(population,"simulated",verbose=verbose,debugMode=debugMode)
    sink()
    file.remove(aa)
    e1 <- proc.time()
    cross <- rearrangeMarkers(cross, population, populationType, cur_map, threshold, placeUsing, flagged, env, addMarkers=TRUE, keep.redundant, chr, verbose=verbose)
  }
  e1 <- proc.time()
  if(verbose && debugMode==2)cat("Enrichment of original map done in:",(e1-s1)[3],"seconds.\n")
  
  #*******ORDERING NEW MAP*******
  if(use.orderMarkers){
    if(verbose)cat("Ordering markers inside the cross object\n")
    s1 <- proc.time()
    aa <- tempfile()
    sink(aa)
    cross <- orderMarkers(cross,use.ripple=FALSE,verbose=TRUE)
    sink()
    file.remove(aa)
    e1 <- proc.time()
    if(verbose && debugMode==2)cat("Saving data into cross object done in:",(e1-s1)[3],"seconds.\n")
   }
   n.newM <- sum(nmar(cross))-  n.originalM
   percentageSat <-  (n.newM/n.originalM)*100
   if(verbose) cat("\ncross.saturate statistics:\n # original markers:",n.originalM,"\n # inserted markers: ",n.newM,"\n saturation (% of markers added): ",percentageSat,"\n")
  invisible(cross)
}


###########################################################################################################
#                                    *** rearrangeMarkers ***
#
# DESCRIPTION:
#   ordering chromosomes using genetic/physical map and corelation rule
# OUTPUT:
#  object of class cross
############################################################################################################
rearrangeMarkers <- function(cross, population, populationType, cur_map, threshold=3, placeUsing, flagged, env, addMarkers=FALSE, keep.redundant=FALSE, chr, verbose=FALSE){
  if(verbose) cat("old map contains",max(cur_map[,1]),"chromosomes\n")
  redundant <- NULL
  if(placeUsing=="qtl"){
    markersNewPostions <- bestQTL.internal(cross,population,threshold,flagged,verbose)
  }else{
    markersNewPostions <- bestCorelated.internal(cross,population,cur_map,threshold,verbose)
  }
  if(verbose) cat("selected",nrow(markersNewPostions),"markers for further analysis\n")

  returncross <- cross
  returncross$geno <- vector(length(unique(cur_map[,1])), mode="list")
  returncross$pheno <- pull.pheno(cross)
  oldnames_ <- rownames(cur_map)[which(rownames(cur_map) %in% rownames(population$offspring$genotypes$real))]
  if(verbose) cat("Reordering markers \n")
  left <- 0
  for(x in 1:length(returncross$geno)){
    #if(verbose) cat("- chr ",x," -\n")    
    oldnamesChr <- rownames(cur_map)[which(cur_map[,1]==x)]
    oldnames <- oldnamesChr[which(oldnamesChr %in% oldnames_)]
    oldpositions <- cur_map[oldnames,2]
    if(x %in% chr){
      newnames_ <- rownames(markersNewPostions)[which(markersNewPostions[,1]==x)]             #TODO: Never Never Never Never Identifiers with _
      if(any(newnames_%in%oldnames)) newnames_ <- newnames_[-which(newnames_%in%oldnames)]
      newnames <- NULL
      newpositions <- NULL
        positions <- cbind(newnames_,markersNewPostions[newnames_,2],markersNewPostions[newnames_,3])
        for(pos in unique(positions[,2])){
          mappingM <- which(positions[,2]==pos)
          newpositions <- c(newpositions,pos)
          if(length(mappingM)==1){
            newnames <- c(newnames,positions[mappingM,1])
          }else{
            selM <- positions[mappingM,1]
            bestM <- which.max(as.numeric(positions[selM,3]))
            left <- left+1
            newnames <- c(newnames,selM[bestM])
          }
        }
    }else{
      newnames <- NULL
      newpositions <- NULL
    }
    toRmv <- NULL
    if(length(newnames)>0 && length(toRmv)>0){
      newnames <- newnames[-toRmv]
      newpositions <- newpositions[-toRmv]
      redundant <- c(redundant,newnames[toRmv])
    }
    if(x %in% chr) if(verbose) cat("Selected:",length(newnames),"new and",length(oldnames),"original markers,",length(toRmv),"markers were removed\n") 
    if(addMarkers){
      returncross$geno[[x]]$data <- insertMarkers.internal(pull.geno(cross)[,newnames],newpositions,t(population$offspring$genotypes$real[oldnames,]),oldpositions, env, populationType)
      newmap <- c(as.numeric(newpositions),oldpositions)
      names(newmap) <- c(newnames,oldnames)
      newmap <- sort(newmap)
      colnames(returncross$geno[[x]]$data) <- c(newnames,oldnames)
      returncross$geno[[x]]$data <- returncross$geno[[x]]$data[,names(newmap)]
    }else{
      returncross$geno[[x]]$data <- pull.geno(cross)[,newnames]
      newmap <- as.numeric(newpositions)
      names(newmap) <- newnames
      newmap <- sort(newmap)
      colnames(returncross$geno[[x]]$data) <- newnames
      returncross$geno[[x]]$data <- returncross$geno[[x]]$data[,names(newmap)]
    }
    returncross$geno[[x]]$map <- c(newmap)
  }
  names(returncross$geno) <- 1:length(returncross$geno)
  for(i in 1:length(returncross$geno)){
    class(returncross$geno[[i]]) <- "A"
  }
  returncross$redundant <- redundant
  returncross$left <- left
  invisible(returncross)
}

###
insertMarkers.internal <- function(newgeno,newpositions,oldgeno,oldpositions,env,populationType){
  if(length(newgeno)<1){ return(oldgeno) }
  toRmv <- NULL    #TODO: Give a description of what I do
  toInv <- NULL    #TODO: Give a description of what I do
  if(is.null(dim(newgeno))){ newgeno <- as.matrix(newgeno) }  #Does this do anything ? If there is no dim how would as.matrix figure it out then ?
  if(is.null(dim(oldgeno))){ oldgeno <- as.matrix(oldgeno) }  #Does this do anything ? If there is no dim how would as.matrix figure it out then ?

  for(i in 1:length(newpositions)){
    distance <- abs(oldpositions-as.numeric(newpositions[i]))
    curCor <- cor(newgeno[,i],oldgeno[,which.min(distance)],use="pair")
    if(abs(curCor) < 0.1){
      toRmv <- c(toRmv,i)
    }else if(curCor < (-0.3)){
      toInv <- c(toInv,i)
    }else{
      #TODO: Figure out what we need to do else
    }
  }
  
  if(populationType == "f2"){ #TODO: Updated this very primitive inversion
    invertM <- newgeno[,toInv]
    invertM[which(invertM==1)] <- 3
    invertM[which(invertM==3)] <- 1
    invertM[which(invertM==5)] <- 4
    invertM[which(invertM==4)] <- 5
    newgeno[,toInv] <- invertM
  }else{ newgeno[,toInv] <- 3 - newgeno[,toInv] }
  
  return(cbind(newgeno,oldgeno))
}

cleanGeno.internal <- function(genoCol,env,genos){
  for(envVal in unique(env)){
    incorr <- 0
    for(geno in genos){
      if(sum(which(genoCol==geno) %in% which(env==envVal)) < 4) incorr <- 1
    }
    if(incorr!=0) genoCol[which(env==envVal)] <- NA
  }
  invisible(genoCol)
}

###########################################################################################################
#                                    *** bestCorelated.internal ***
#
# DESCRIPTION:
#   subfunction of segragateChromosomes.internal, returns matrix showing for every reco map chromosome from 
#  which physicall map chromosome majority of markers comes
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
############################################################################################################
bestCorelated.internal <- function(cross, population, cur_map, corSDTreshold, verbose=FALSE){
  cormatrix <- map2mapCorrelationMatrix(cross,population,verbose)
  maximums <- apply(abs(cormatrix), 2,max)
  means <- apply(abs(cormatrix), 2,mean)
  sds <- apply(abs(cormatrix), 2,sd)
  selected <- which(maximums > (means+corSDTreshold*sds))  # Select markers that are correlated highly with more than one of the old markers
  cormatrix <- cormatrix[,selected]
  bestCorMarkers <- matrix(0,length(selected),2)
  bestCorMarkers[,1] <- apply(abs(cormatrix),2,function(r){rownames(cormatrix)[which.max(r)]})
  bestCorMarkers[,2] <- apply(abs(cormatrix),2,function(r){rownames(cormatrix)[which.max(r[-which.max(r)])]})
  bestCorMarkers[,3] <- apply(abs(cormatrix),2,max)
  rownames(bestCorMarkers) <- rownames(cormatrix)
  output <- t(apply(bestCorMarkers,1,bestCorelatedSub.internal,cur_map))
  invisible(output)
}

bestCorelatedSub.internal <- function(bestCorMarkersRow,cur_map){
  chr <- cur_map[bestCorMarkersRow[1],1]
  pos <- mean(cur_map[bestCorMarkersRow[1],2],cur_map[bestCorMarkersRow[2],2])
  invisible(c(chr,pos,bestCorMarkersRow[3]))
}

bestQTLSub.internal <- function(qtls,marker){
  cur_max <- which.max(qtls$lod[marker,])
  cur_row <- c(qtls$chr[marker,cur_max], qtls$pos[marker,cur_max], max(qtls$lod[marker,]))
  invisible(cur_row)
}

###########################################################################################################
#                                    *** bestQTL.internal ***
#
# DESCRIPTION:
#   subfunction of segragateChromosomes.internal, returns matrix showing for every reco map chromosome from 
#  which physicall map chromosome majority of markers comes
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
############################################################################################################
bestQTL.internal <- function(cross, population, threshold, flagged, verbose=FALSE){
  genotypes <- population$offspring$genotypes$real
  if(is.null(population$offspring$genotypes$qtl$flags)) stop("Old version of the QTL scan detected. Re-run scan.qtls!")
  markers <- markernames(cross)
  phenotypes <- pull.geno(cross)[,markers]
  output <- NULL
  count <- 0
  if(verbose) cat("Starting qtl analysis.\n")
  s<- proc.time()
  if(is.null(population$offspring$genotypes$qtl)) stop("No qtl data in population$offspring$genotypes$qtl, run scan.qtls function first.")
  peaksMatrix <- getpeaks.internal(abs(population$offspring$genotypes$qtl$lod),threshold)
  e<- proc.time()
  envInf <- 0
  epiInf <- 0
  if(verbose) cat("Qtl analysis done in:",(e-s)[3],"seconds\n")
  rownames(peaksMatrix) <- markers
  for(marker in markers){
    if(sum(peaksMatrix[marker,]==2)==1){ #TODO: Figure out the logic here, Its not logical
      if(any(population$offspring$genotypes$qtl$flags[marker,] > (threshold/2))){
        envInf <- envInf + 1
        if(flagged=="remove"){
          cat("Marker:",marker,"shows significant association with environent and will be removed.\n")
          output <- rbind(output,c(NA,NA,NA))
        }else{
          if(flagged=="warn") cat("Marker:",marker,"shows significant association with environent.\n")
          output <- rbind(output,bestQTLSub.internal(population$offspring$genotypes$qtl,marker))
        }
      }else if(population$offspring$genotypes$qtl$logLik[marker,1] > (population$offspring$genotypes$qtl$logLik[marker,2]+10)){
        epiInf <- epiInf + 1
        if(flagged=="remove"){
          cat("Marker:",marker,"is influenced by an epistatic interaction and will be removed.\n")
          output <- rbind(output,c(NA,NA,NA))
        }else{
          if(flagged=="warn") cat("Marker:",marker,"is influenced by an epistatic interaction.\n")
          output <- rbind(output,bestQTLSub.internal(population$offspring$genotypes$qtl,marker))
        }
      }else{
        output <- rbind(output,bestQTLSub.internal(population$offspring$genotypes$qtl,marker))
      }
    }else{
      count <- count+1
      output <- rbind(output,c(NA,NA,NA))
    }
  }
  #to have same format of the output as in bestcorrelated
  if(verbose && flagged=="remove"){
    cat("Removed:",envInf,"markers showing significant association with environent.\n")
    cat("Removed:",epiInf,"markers influenced by an epistatic interaction.\n")
  }else if(verbose){
    cat(envInf,"markers show significant association with environent.\n")
    cat(epiInf,"markers are influenced by an epistatic interaction.\n")
  }
  rownames(output) <- markers
  if(any(is.na(output[,1]))) output <- output[-which(is.na(output[,1])),]
  if(any(is.na(output[,2]))) output <- output[-which(is.na(output[,2])),]
  invisible(output)
}

#TODO: Add documentation
fullScanRow.internal <- function(genoRow, phenoRow, env){
  model <- lm(phenoRow ~ env + genoRow + env:genoRow)
  return(c(-log10(anova(model)[[5]])[1:3], logLik(model)))
}

###########################################################################################################
#                                    *** QTLscan.internal ***
#
# DESCRIPTION:
# subfunction by Danny Arends to map QLTs modfied to work on a single phenotype
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
############################################################################################################
scan.qtls <- function(population,map=c("genetic","physical"), env, step=0.1,verbose=FALSE){
  if(missing(population)) stop("Please provide a population object\n")
  check.population(population)
  if(is.null(population$offspring$genotypes$real)){
    stop("No original genotypes in population$offspring$genotypes$real, load them in using add.to.population\n")
  }
  if(is.null(population$offspring$genotypes$simulated)){
    stop("No simulated genotypes in population$offspring$genotypes$simulated, run generate.biomarkers first\n")
  }
  if(missing(env)) env <- rep(1,ncol(population$offspring$phenotypes))
  map <- match.arg(map)
  if(map=="genetic"){
    matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$genetic))
    if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
    if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
      population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
      population$maps$genetic <- population$maps$genetic[rownames(population$offspring$genotypes$real),]
      n.markersToRmv <- nrow(population$offspring$genotypes$real)-length(matchingMarkers)
      if(verbose && n.markersToRmv>0) cat(n.markersToRmv,"markers were removed due to name mismatch\n")
    }
    population10pheno <- population
    population10pheno$offspring$phenotypes <- population10pheno$offspring$phenotypes[1:10,]
    aa <- tempfile()
    sink(aa)          #TODO: When using Sink make sure you Try{}Catch everything, we need to dis-able sink even if everythign exploded
    returncross <- genotypesToCross.internal(population10pheno,"real","map_genetic")
    sink()
    file.remove(aa)   #TODO: If we have an error we don't delete our file ????
  }else{
    matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$physical))
    if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
    if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
      population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
      population$maps$physical <- population$maps$physical[rownames(population$offspring$genotypes$real),]
      n.markersToRmv <- nrow(population$offspring$genotypes$real)-length(matchingMarkers)
      if(verbose && n.markersToRmv>0) cat(n.markersToRmv,"markers were removed due to name mismatch\n")
    }
    #TODO: Why is this here the original comment: 'for faster creation of cross' is meaningless
    population10pheno <- population
    population10pheno$offspring$phenotypes <- population10pheno$offspring$phenotypes[1:10,]
    aa <- tempfile()
    sink(aa)
    returncross <- genotypesToCross.internal(population10pheno,"real","map_physical")
    sink()
    file.remove(aa)
  }
  returncross$pheno <- t(population$offspring$genotypes$simulated)
  returncross <- calc.genoprob(returncross,step=step)
  returncrosstwo <- calc.genoprob(returncross,step=0)
  if(verbose) cat("Starting qtl scan, this may take a long time to finish!\n")
  s <- proc.time()
  lod <- NULL
  pos <- NULL
  chr <- NULL
  flags <- NULL
  names_ <- NULL
  scan2 <- NULL
  logLikeli <- NULL
  done <- 0
  useEnv <- TRUE
  if(!(length(unique(env))>1)) useEnv <- FALSE
  for(i in 1:nrow(population$offspring$genotypes$simulated)){
    curlogLikeli <- NULL
    phenotype <- pull.pheno(returncross)[,i]
    perc <- round(i*100/nrow(population$offspring$genotypes$simulated))
    if(perc%%10==0 && !(perc%in%done)){
      e <- proc.time()
      cat("Analysing markers",perc,"% done, estimated time remaining:",(e-s)[3]/perc*(100-perc),"s\n")
      done <- c(done,perc)
    }
    if(useEnv){
      flag <- t(apply(pull.geno(returncross),2,fullScanRow.internal,phenotype,env))
      flags <- rbind(flags,c(max(flag[,1]),max(flag[,3])))
      curlogLikeli <- c(curlogLikeli,min(flag[,4]))
    }else{
      flags <- rbind(flags,c(0,0))
      curlogLikeli <- c(curlogLikeli,0)
    }
    curScan <- scanone(returncross,pheno.col=i,model="np")
    aa <- tempfile()                #TODO:  When using Sink make sure you Try{}Catch everything, we need to dis-able sink even if everythign exploded
    sink(aa)
    curScantwo <- scantwo(returncrosstwo,pheno.col=i)
    maxLine <- which.max(summary(curScantwo)[,6])
    chr1 <- summary(curScantwo)[maxLine,1]
    marker1 <- which(returncross$geno[[chr1]]$map==summary(curScantwo)[maxLine,3])
    chr2 <- summary(curScantwo)[maxLine,2]
    marker2 <- which(returncross$geno[[chr2]]$map==summary(curScantwo)[maxLine,4])
    genoRow1 <- returncross$geno[[chr1]]$data[,marker1]
    genoRow2 <- returncross$geno[[chr2]]$data[,marker2]
    model <- lm(phenotype ~ genoRow1 + genoRow2 + genoRow1:genoRow2)
    curlogLikeli <- c(curlogLikeli,logLik(model))
    sink()
    file.remove(aa)

    chr <- rbind(chr,curScan[,1])
    pos <- rbind(pos,curScan[,2])
    lod <- rbind(lod,curScan[,3])
    logLikeli <- rbind(logLikeli,curlogLikeli)
    names_ <- c(names_,colnames(returncross$pheno)[i])
  }

  e <- proc.time()
  if(verbose) cat("Qtl scan done in",(e-s)[3],"s\n")
  population$offspring$genotypes$qtl$lod <- lod
  population$offspring$genotypes$qtl$pos <- pos
  population$offspring$genotypes$qtl$chr <- chr
  population$offspring$genotypes$qtl$flags <- flags
  population$offspring$genotypes$qtl$logLik <- logLikeli
  population$offspring$genotypes$qtl$names <- names_
  rownames(population$offspring$genotypes$qtl$lod) <- names_
  colnames(population$offspring$genotypes$qtl$lod) <- rownames(curScan)   
  rownames(population$offspring$genotypes$qtl$chr) <- names_
  colnames(population$offspring$genotypes$qtl$chr) <- rownames(curScan)  
  rownames(population$offspring$genotypes$qtl$pos) <- names_
  colnames(population$offspring$genotypes$qtl$pos) <- rownames(curScan)
  rownames(population$offspring$genotypes$qtl$logLik) <- names_
  rownames(population$offspring$genotypes$qtl$flags) <- names_
  invisible(population)
}

###########################################################################################################
#                                    *** QTLscan.internal ***
#
# DESCRIPTION:
# subfunction by Danny Arends to map QLTs modfied to work on a single phenotype
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
############################################################################################################
getpeaks.internal <- function(qtlprofiles, cutoff = 4.0){
  if(!any(qtlprofiles==Inf)) qtlprofiles[which(qtlprofiles==Inf)] <- 1000
  mmatrix <- NULL
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
    for(a in which(maximums>0)){ mrow[maximums[a]] <- 2 }
    for(b in which(maximums<0)){ mrow[(-maximums[b])] <- -2 }
    mmatrix <- rbind(mmatrix,mrow)
  }
  mmatrix
}

###########################################################################################################
#                                    *** map2mapCorrelationMatrix.internal ***
#
# DESCRIPTION:
#   calculating correlation matrix between genotypes inside cross object and ones from population
# OUTPUT:
#   matrix of correlations
############################################################################################################
map2mapCorrelationMatrix<- function(cross,population,verbose=FALSE){
  if(missing(cross)) stop("Please provide a cross object\n")
  if(missing(population)) stop("Please provide original genotypes\n")

  genotypes <- pull.geno(cross)
  if(verbose) cat("Calculating correlation matrix\n")
  if(!is.null(population$offspring$genotypes$real)){
    genotypesCorelationMatrix <- apply(genotypes,2,function(cgc){cor(cgc,t(population$offspring$genotypes$real),use="pair")})
    colnames(genotypesCorelationMatrix) <- colnames(genotypes)
    rownames(genotypesCorelationMatrix) <- rownames(population$offspring$genotypes$real)    
    invisible(genotypesCorelationMatrix)
  }else{
    stop("Load known genotypes into the population using add.to.population(p,genotypes,\"offspring$genotypes\")")
  }
}

###########################################################################################################
#                                           *** map2mapImage ***
#
# DESCRIPTION:
#   subfunction of segragateChromosomes.internal, returns matrix showing for every reco map chromosome from 
#  which physicall map chromosome majority of markers comes# 
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
############################################################################################################
map2mapImage <- function(genotypesCorelationMatrix,population,cross,corThreshold=0.5,verbose=FALSE){
  if(missing(genotypesCorelationMatrix)){
    cat("Correlation matrix not provided, calulating one")
    genotypesCorelationMatrix <- map2mapCorrelationMatrix(cross,population,verbose)
    if(missing(cross)) stop("No object of class cross, please run either cross.denovo or enrichExistingMap\n")
    if(missing(population)) stop("Please provide a population object\n")
    check.population(population)
  }
  heatmap(genotypesCorelationMatrix,breaks = c(-1,-corThreshold,corThreshold,1),col=c("blue","white","red"),Rowv=NA)
}

