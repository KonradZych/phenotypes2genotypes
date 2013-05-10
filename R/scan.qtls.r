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
  
  if(is.null(population$offspring$genotypes$real))      stop("No original genotypes in population$offspring$genotypes$real, load them in using add.to.population\n")
  if(is.null(population$offspring$genotypes$simulated)) stop("No simulated genotypes in population$offspring$genotypes$simulated, run generate.biomarkers first\n")
  
  if(missing(env)) env <- rep(1,ncol(population$offspring$phenotypes)) #if there is no infor about env -> all of them in the same env
  
  map <- match.arg(map)

  if(map=="genetic"){
    population      <- matchMarkers(population, population$maps$genetic, mapType="genetic")
    originalMap     <- population$maps$genetic
  }else{
    population      <- matchMarkers(population, population$maps$physical, mapType="physical")
    originalMap     <- population$maps$physical
  }
  
  population10pheno                      <- population
  population10pheno$offspring$phenotypes <- population10pheno$offspring$phenotypes[1:10,]
  
  ### creation of the cross so that we can use r/qtl for qtl mapping
  tryCatch({
    aa <- tempfile()
    sink(aa)
    returncross <- genotypesToCross.internal(population10pheno,"real","map_genetic")
  },
  error= function(err){
    print(paste("ERROR in scan.qtls while creating cross:  ",err))
    sink()            # sink if errored -> otherwise everything is sinked into aa file
    # file is not removed -> contains output that may help with debugging
  },
  finally={
    sink()
    file.remove(aa) # no error -> close sink and remove unneeded file
  })
  
  returncross$pheno <- t(population$offspring$genotypes$simulated)
  returncross       <- calc.genoprob(returncross,step=step)
  s                 <- proc.time()
  lod               <- NULL # LOD scores
  pos               <- NULL # position of the QTL peak on a chromosome
  chr               <- NULL # chromosome of the peak marker
  interactions      <- NULL # 1 - max lod score for Env, 2 - G:E, 3- epistasis
  markerNames       <- NULL # names of the markers se
  done              <- 0    
  useEnv            <- TRUE

  if(!(length(unique(env))>1)) useEnv <- FALSE

  for(i in 1:nrow(population$offspring$genotypes$simulated)){
    curlogLikeli <- NULL
    phenotype <- pull.pheno(returncross)[,i]
    perc <- round(i*100/nrow(population$offspring$genotypes$simulated))
    if(perc%%10==0 && !(perc%in%done)){ # bit dirty hack to avoid displaying the same value more than once in verbose mode (this can happen due to rounding)
      e <- proc.time()
      cat("Analysing markers",perc,"% done, estimated time remaining:",(e-s)[3]/perc*(100-perc),"s\n")
      done <- c(done,perc)
    }
    if(useEnv){
      fullScanResults    <- t(apply(pull.geno(returncross),2,fullScanRow.internal,phenotype,env))
      curInteractions    <- c(max(fullScanResults[,1]),max(fullScanResults[,3]))
    }else{
      curInteractions    <- c(0,0)
    }
    tryCatch({
      aa <- tempfile()
      sink(aa)
      curScan <- scanone(returncross,pheno.col=i,model="np")
    },
    error= function(err){
      print(paste("ERROR in scan.qtls while using scanone:  ",err))
      sink()            # sink if errored -> otherwise everything is sinked into aa file
      # file is not removed -> contains output that may help with debugging
    },
    finally={
      sink()
      file.remove(aa) # no error -> close sink and remove unneeded file
    })
    
    #epistaticInter  <- checkForEpistasis(curScan,pull.geno(returncross),pull.pheno(returncross)[,i],env)
    curScantwo      <- scantwo(returncrosstwo,pheno.col=i)
    epistaticInter  <- max(summary(curScantwo)[,6])
    curInteractions <- c(curInteractions,epistaticInter)
    interactions    <- rbind(interactions,curInteractions)

    chr             <- rbind(chr,curScan[,1])
    pos             <- rbind(pos,curScan[,2])
    lod             <- rbind(lod,curScan[,3])
    markerNames     <- c(markerNames,colnames(returncross$pheno)[i])
  }

  e <- proc.time()
  if(verbose) cat("Qtl scan done in",(e-s)[3],"s\n")
  
  ### packing all the results into the population object
  population$offspring$genotypes$qtl$lod                     <- lod
  rownames(population$offspring$genotypes$qtl$lod)           <- markerNames
  colnames(population$offspring$genotypes$qtl$lod)           <- rownames(curScan)
  
  population$offspring$genotypes$qtl$pos                     <- pos
  rownames(population$offspring$genotypes$qtl$chr)           <- markerNames
  colnames(population$offspring$genotypes$qtl$chr)           <- rownames(curScan)

  population$offspring$genotypes$qtl$chr                     <- chr
  rownames(population$offspring$genotypes$qtl$pos)           <- markerNames
  colnames(population$offspring$genotypes$qtl$pos)           <- rownames(curScan)
  
  population$offspring$genotypes$qtl$interactions            <- interactions
  rownames(population$offspring$genotypes$qtl$interactions)  <- markerNames
  invisible(population)
}

checkForEpistasis <- function(scanResults,originalGeno,marker,env){
  maxGenoNr <- rownames(scanResults)[which.max(scanResults[,3])]
  results   <- apply(originalGeno, 2, twoGenosModel, marker, originalGeno[,maxGenoNr], env)
  results   <- results[-maxGenoNr] # removing model with env + maxGenoNr +maxGenoNr (overfit!)
  invisible(max(results))
}

twoGenosModel <- function(genoRow, markerRow, maxGeno, env){
  curRes <- -log10(anova(lm(markerRow~env+maxGeno+genoRow))[[5]][3])
  invisible(curRes)
}

#TODO: Add documentation
fullScanRow.internal <- function(genoRow, phenoRow, env){
  model <- lm(phenoRow ~ env + genoRow + env:genoRow)
  return(c(-log10(anova(model)[[5]])[1:3], logLik(model)))
}

