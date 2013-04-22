#
# scan.qtls.R
#
# Copyright (c) 2010-2013 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified April, 2013
# first written Mar, 2011
# Contains: scan.qtls
#

##########################################################################################################
#                                    *** scan.qtls ***
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
  population$offspring$genotypes$qtl$lod              <- lod
  population$offspring$genotypes$qtl$pos              <- pos
  population$offspring$genotypes$qtl$chr              <- chr
  population$offspring$genotypes$qtl$flags            <- flags
  population$offspring$genotypes$qtl$logLik           <- logLikeli
  population$offspring$genotypes$qtl$names            <- names_
  rownames(population$offspring$genotypes$qtl$lod)    <- names_
  colnames(population$offspring$genotypes$qtl$lod)    <- rownames(curScan)   
  rownames(population$offspring$genotypes$qtl$chr)    <- names_
  colnames(population$offspring$genotypes$qtl$chr)    <- rownames(curScan)  
  rownames(population$offspring$genotypes$qtl$pos)    <- names_
  colnames(population$offspring$genotypes$qtl$pos)    <- rownames(curScan)
  rownames(population$offspring$genotypes$qtl$logLik) <- names_
  rownames(population$offspring$genotypes$qtl$flags)  <- names_
  invisible(population)
}

#TODO: Add documentation
fullScanRow.internal <- function(genoRow, phenoRow, env){
  model <- lm(phenoRow ~ env + genoRow + env:genoRow)
  return(c(-log10(anova(model)[[5]])[1:3], logLik(model)))
}
