#
# cross.saturate.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
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
cross.saturate <- function(population, cross, map=c("genetic","physical"), placeUsing=c("qtl","correlation"), threshold=3, chr, use.orderMarkers=FALSE, gffFile, verbose=FALSE, debugMode=0){
  if(missing(population)) stop("Please provide a population object\n")
  populationType <- class(population)[2]
  check.population(population)
  if(!is.numeric(threshold)||is.na(threshold)) stop("Please provide correct threshold")
  if(threshold<=0){
    cat("WARNING: threshold too low, all possible markers will be selected\n")
  }else if(placeUsing=="correlation" && threshold>=5){
    cat("WARNING: threshold too high, few new markers will be selected\n")
  }else if(placeUsing=="qtl" && threshold>=20){
    cat("WARNING: threshold too high, few new markers will be selected\n")
  }
  map <- match.arg(map)
  placeUsing <- checkParameters.internal(placeUsing,c("qtl","correlation"),"placeUsing")
  if(missing(cross)){
    if(is.null(population$offspring$genotypes$real)){
      stop("No original genotypes in population$offspring$genotypes$real, load them in using add.to.population\n")
    }
    if(is.null(population$offspring$genotypes$simulated)){
      stop("No genotype data in population$offspring$genotypes$simulated, run generate.biomarkers first\n")
    }else{
      #*******SAVING CROSS OBJECT*******
      s1 <- proc.time()
      aa <- tempfile()
      sink(aa)
      cross <- genotypesToCross.internal(population,"simulated",verbose=verbose,debugMode=debugMode)
      sink()
      file.remove(aa)
      e1 <- proc.time()
      if(verbose && debugMode==2)cat("Saving data into cross object done in:",(e1-s1)[3],"seconds.\n")
    }
  }else{
    population <- set.geno.from.cross(cross,population,map)
    population <- scan.qtls(population,map)
      aa <- tempfile()
      sink(aa)
      cross <- genotypesToCross.internal(population,"simulated",verbose=verbose,debugMode=debugMode)
      sink()
      file.remove(aa)
      e1 <- proc.time()
  }
  if(!(all(rownames(population$offspring$genotypes$simulated)%in%rownames(population$offspring$genotypes$qtl$lod)))){
    stop("QTL scan results don't match with simulated genotypes, please, run scan.qtls function\n")
  }else if(!(all(rownames(population$offspring$genotypes$qtl$lod)%in%rownames(population$offspring$genotypes$simulated)))){
    stop("QTL scan results don't match with simulated genotypes, please, run scan.qtls function\n")
  }
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
  cross <- rearrangeMarkers(cross, population, populationType, cur_map, threshold, placeUsing, addMarkers=TRUE, chr, gffFile, verbose=verbose)
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
rearrangeMarkers <- function(cross, population, populationType, cur_map, threshold=3, placeUsing, addMarkers=FALSE, chr, gffFile, verbose=FALSE){
  if(verbose) cat("old map contains",max(cur_map[,1]),"chromosomes\n")
  nonReduntant <- NULL
  if(placeUsing=="qtl"){
    markersNewPostions <- bestQTL.internal(cross,population,threshold,verbose)
  }else{
    markersNewPostions <- bestCorelated.internal(cross,population,cur_map,threshold,verbose)
  }
  if(nrow(markersNewPostions) == 0){
    cat("selected",nrow(markersNewPostions),"markers with current corThreshold, there will be only markers from old map in the cross object\n")
  }else if(verbose){
    cat("selected",nrow(markersNewPostions),"markers for further analysis\n")
  }
  returncross <- cross
  returncross$geno <- vector(length(unique(cur_map[,1])), mode="list")
  returncross$pheno <- pull.pheno(cross)
  oldnames_ <- rownames(cur_map)[which(rownames(cur_map) %in% rownames(population$offspring$genotypes$real))]
  if(verbose) cat("Reordering markers \n")  
  for(x in 1:length(returncross$geno)){
    #if(verbose) cat("- chr ",x," -\n")    
    oldnamesChr <- rownames(cur_map)[which(cur_map[,1]==x)]
    oldnames <- oldnamesChr[which(oldnamesChr %in% oldnames_)]
    oldpositions <- cur_map[oldnames,2]
    if(x %in% chr){
    newnames <- rownames(markersNewPostions)[which(markersNewPostions[,1]==x)]
    if(any(newnames%in%oldnames)){
      newnames <- newnames[-which(newnames%in%oldnames)]
    }
    newpositions <- markersNewPostions[newnames,2]
  }else{
    newnames <- NULL
    newpositions <- NULL
  }
    toRmv <- NULL
    if(length(newnames)>0){
      for(i in 1:length(newpositions)){
        if(newpositions[i]%in%oldpositions){
          toRmv <- c(toRmv,i)
        }
      }
      if(length(toRmv)>0){
        newnames <- newnames[-toRmv]
        newpositions <- newpositions[-toRmv]
        nonReduntant <- c(nonReduntant,newnames)
        }
    }
     if(x %in% chr) if(verbose) cat("Selected:",length(newnames),"new and",length(oldnames),"original markers,",length(toRmv),"markers were removed\n") 
    if(addMarkers){
      returncross$geno[[x]]$data <- insertMarkers.internal(pull.geno(cross)[,newnames],newpositions,t(population$offspring$genotypes$real[oldnames,]),oldpositions,  populationType)
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
  if(!(missing(gffFile))){
    cat("Saving gff files.\n")
    filename1 <- paste(gffFile,"all.gff",sep="")
    filename2 <- paste(gffFile,"nonRedundant.gff",sep="")
      newnames <- rownames(markersNewPostions)
      if(!(any(newnames %in% rownames(population$maps$physical)))){
        if(!(any(oldnames_%in% rownames(population$maps$physical)))){
          cat("population object doesn't contain information about positions of the markers, gff file won't be saved\n")
        }else{
          markers <- population$offspring$genotypes$real[oldnames_,]
          positions <- population$maps$physical[rownames(markers),]
          saveGff.internal(gffFile,markers,positions)
        }
      }else if(!(any(oldnames_%in% rownames(population$maps$physical)))){
        markers <- t(pull.geno(cross)[,newnames])
        positions <- population$maps$physical[rownames(markers),]
        saveGff.internal(gffFile,markers,positions,nonReduntant)
        markers <- t(pull.geno(cross)[,nonReduntant])
        positions <- population$maps$physical[rownames(markers),]
        saveGff.internal(filename2,markers,positions)
      }else{
        chrL <- chromosomesLengths.internal(population$maps$physical)
        markers <- rbind(t(pull.geno(cross)[,newnames]),population$offspring$genotypes$real[oldnames_,])
        positions <- population$maps$physical[rownames(markers),]
        saveGff.internal(gffFile,markers,positions,nonReduntant)
        markers <- rbind(t(pull.geno(cross)[,nonReduntant]),population$offspring$genotypes$real[oldnames_,])
        positions <- population$maps$physical[rownames(markers),]
        saveGff.internal(filename2,markers,positions)
      }
  }
  invisible(returncross)
}

###
saveGff.internal <- function(gffFile="population.gff", markers, positions, nonReduntant){
 print(dim(markers))
  cat("##gff-version 3\n",file=gffFile,append=FALSE)
  for(marker in 1:nrow(markers)-1){
    if(marker %in% nonReduntant){cat("Chr",positions[marker,1],"\t.\tmarker\t",positions[marker,2],"\t",positions[marker,3],"\t100\t+\t.\tID=recombination\n",file=gffFile,append=TRUE,sep='')
    }else{cat("Chr",positions[marker,1],"\t.\tredundantmarker\t",positions[marker,2],"\t",positions[marker,3],"\t100\t+\t.\tID=recombination\n",file=gffFile,append=TRUE,sep='')}
    #if(marker %in% nonReduntant)cat("Chr",positions[marker,1],"\t.\tnonreduntantmarker\t",(positions[marker,2]+positions[marker+1,2])/2,"\t",(positions[marker,2]+positions[marker+1,2])/2,"\t100\t+\t.\tID=recombination\n",file=gffFile,append=TRUE,sep='')
  }
}

insertMarkers.internal <- function(newgeno,newpositions,oldgeno,oldpositions,populationType){
  if(length(newgeno)<1){
    return(oldgeno)
  }
  toRmv <- NULL
  toInv <- NULL
  if(is.null(dim(newgeno))){
    newgeno <- as.matrix(newgeno)
  }
  if(is.null(dim(oldgeno))){
    oldgeno <- as.matrix(oldgeno)
  }
  for(i in 1:length(newpositions)){
    distance <- abs(oldpositions-newpositions[i])
    curCor <- cor(newgeno[,i],oldgeno[,which.min(distance)],use="pair")
    #print(curCor)
    #print(toInv)
    #cat(i,":",which.min(distance),":",curCor,"\n")
    if(curCor<0.5 && curCor>(-0.3)){
      toRmv <- c(toRmv,i)
    }else if(curCor<(-0.3)){
      toInv <- c(toInv,i)
    }
  }
  #if(!is.null(toRmv)){
  #  newgeno <- newgeno[,-toRmv]
  #}
  #print(toInv)
  ### very primitive inversion in here!
  if(populationType=="f2"){
    invertM <- newgeno[,toInv]
    invertM[which(invertM==1)] <- 3
    invertM[which(invertM==3)] <- 1
    invertM[which(invertM==5)] <- 4
    invertM[which(invertM==4)] <- 5
    newgeno[,toInv] <- invertM
  }else{
    newgeno[,toInv] <- 3 - newgeno[,toInv]
  }
  return(cbind(newgeno,oldgeno))
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
bestCorelated.internal <- function(cross,population, cur_map,corSDTreshold,verbose=FALSE){
  cormatrix <- map2mapCorrelationMatrix(cross,population,verbose)
  maximums <- apply(abs(cormatrix),2,max)
  means <- apply(abs(cormatrix),2,mean)
  sds <- apply(abs(cormatrix),2,sd)
  #select markers that are correlated highly with more than one of the old markers
  selected <- which(maximums > (means+corSDTreshold*sds))
  cormatrix <- cormatrix[,selected]
  bestCorMarkers <- matrix(0,length(selected),2)
  bestCorMarkers[,1] <- apply(abs(cormatrix),2,function(r){rownames(cormatrix)[which.max(r)]})
  bestCorMarkers[,2] <- apply(abs(cormatrix),2,function(r){rownames(cormatrix)[which.max(r[-which.max(r)])]})
  rownames(bestCorMarkers) <- rownames(cormatrix)
  output <- t(apply(bestCorMarkers,1,bestCorelatedSub.internal,cur_map))
  invisible(output)
}

bestCorelatedSub.internal <- function(bestCorMarkersRow,cur_map){
  chr <- cur_map[bestCorMarkersRow[1],1]
  pos <- mean(cur_map[bestCorMarkersRow[1],2],cur_map[bestCorMarkersRow[2],2])
  invisible(c(chr,pos))
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
bestQTL.internal <- function(cross, population, treshold,verbose=FALSE){
  genotypes <- population$offspring$genotypes$real
  markers <- markernames(cross)
  phenotypes <- pull.geno(cross)[,markers]
  output <- NULL
  if(verbose) cat("Starting qtl analysis.\n")
  s<- proc.time()
  if(is.null(population$offspring$genotypes$qtl)) stop("No qtl data in population$offspring$genotypes$qtl, run scan.qtls function first.")
  peaksMatrix <- getpeaks.internal(abs(population$offspring$genotypes$qtl$lod),treshold)
  e<- proc.time()
  if(verbose) cat("Qtl analysis done in:",(e-s)[3],"seconds\n")
  rownames(peaksMatrix) <- markers
  for(marker in markers){
    if(sum(peaksMatrix[marker,]==2)==1){
      cur_max <- which.max(population$offspring$genotypes$qtl$lod[marker,])
      cur_row <- c(population$offspring$genotypes$qtl$chr[marker,cur_max],population$offspring$genotypes$qtl$pos[marker,cur_max])
      output <- rbind(output,cur_row)
    }else{
      output <- rbind(output,c(NA,NA))
    }
  }
  #to have same format of the output as in bestcorrelated
  rownames(output) <- markers
  if(any(is.na(output[,1]))) output <- output[-which(is.na(output[,1])),]
  if(any(is.na(output[,2]))) output <- output[-which(is.na(output[,2])),]
  invisible(output)
}


###########################################################################################################
#                                    *** QTLscan.internal ***
#
# DESCRIPTION:
# subfunction by Danny Arends to map QLTs modfied to work on a single phenotype
# OUTPUT:
#  vector with new ordering of chromosomes inside cross object
############################################################################################################
scan.qtls <- function(population,map=c("genetic","physical"),step=0.1,verbose=FALSE){
  if(missing(population)) stop("Please provide a population object\n")
  check.population(population)
  if(is.null(population$offspring$genotypes$real)){
    stop("No original genotypes in population$offspring$genotypes$real, load them in using add.to.population\n")
  }
  if(is.null(population$offspring$genotypes$simulated)){
    stop("No simulated genotypes in population$offspring$genotypes$simulated, run generate.biomarkers first\n")
  }
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
    sink(aa)
    returncross <- genotypesToCross.internal(population10pheno,"real","map_genetic")
    sink()
    file.remove(aa)
  }else{
    matchingMarkers <- which(rownames(population$offspring$genotypes$real)%in%rownames(population$maps$physical))
    if(length(matchingMarkers)<=0) stop("Marker names on the map and in the genotypes doesn't match!\n")
    if(length(matchingMarkers)!=nrow(population$offspring$genotypes$real)){
      population$offspring$genotypes$real <- population$offspring$genotypes$real[matchingMarkers,]
      population$maps$physical <- population$maps$physical[rownames(population$offspring$genotypes$real),]
      n.markersToRmv <- nrow(population$offspring$genotypes$real)-length(matchingMarkers)
      if(verbose && n.markersToRmv>0) cat(n.markersToRmv,"markers were removed due to name mismatch\n")
    }
    #for faster creation of cross
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
  if(verbose) cat("Starting qtl scan, this may take a long time to finish!\n")
  s <- proc.time()
  res <- NULL
  pos <- NULL
  chr <- NULL
  names_ <- NULL
  for(i in 1:nrow(population$offspring$genotypes$simulated)){
    if(i%%50==0){
      cat("Analysing marker:",i,"\n")
      }
    curScan <- scanone(returncross,pheno.col=i,method="hk")
    pos <- rbind(pos,curScan[,2])
    res <- rbind(res,curScan[,3])
    chr <- rbind(chr,curScan[,1])
    names_ <- c(names_,colnames(returncross$pheno)[i])
  }
  #population$offspring$genotypes$qtl <- t(matrix(unlist(lapply(markers,QTLscan.internal,phenotypes,genotypes)),nrow(genotypes),length(markers)))
  e <- proc.time()
  if(verbose) cat("Qtl scan done in",(e-s)[3],"s\n")
  population$offspring$genotypes$qtl$lod <- res
  population$offspring$genotypes$qtl$pos <- pos
  population$offspring$genotypes$qtl$chr <- chr
  population$offspring$genotypes$qtl$names <- names_
  rownames(population$offspring$genotypes$qtl$lod) <- names_
  colnames(population$offspring$genotypes$qtl$lod) <- rownames(curScan)   
  rownames(population$offspring$genotypes$qtl$chr) <- names_
  colnames(population$offspring$genotypes$qtl$chr) <- rownames(curScan)  
  rownames(population$offspring$genotypes$qtl$pos) <- names_
  colnames(population$offspring$genotypes$qtl$pos) <- rownames(curScan)
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
  #is.cross(cross)
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
