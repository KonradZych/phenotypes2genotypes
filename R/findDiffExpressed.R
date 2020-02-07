#
# find.diff.expressed.R
#
# Copyright (c) 2010-2012 GBIC: Danny Arends, Konrad Zych and Ritsert C. Jansen
# last modified May, 2012
# first written Mar, 2011
# Contains: find.diff.expressed, fake.population, plotRPpval
#

# find.diff.expressed
#
# DESCRIPTION:
#  Using Rank Product or student t-test analysis to select differentially expressed genes.
# PARAMETERS:
#   - population - Object of class population , must contain founders phenotypic data.
#   - verbose - Be verbose
#   - debugMode - 1: Print our checks, 2: print additional time information
#   - ... - parameters send to RP function
# OUTPUT:
#  An object of class population containing object of class RP in $founders$RP
#
find.diff.expressed <- function(population, verbose=FALSE, debugMode=0, ...){
  #checks
  if(missing(population)) stop("provide population object\n")
  check.population(population)

  if(verbose && debugMode==1) cat("find.diff.expressed starting withour errors in checkpoints.\n")
  
  s<-proc.time()
  population$founders$RP$pval<- t(rbind(apply(population$founders$phenotypes,1,findUsingTTest.internal,population$founders$groups)))
  e<-proc.time()

  if(verbose && debugMode==2)cat("Differentially expressed genes found in:",(e-s)[3],"seconds.\n")
  invisible(population)
}

############################################################################################################
#                  *** findUsingTTest.internal ***
#
# DESCRIPTION:
#  subfunction of find.diff.expressed using t-test to assess whether gene is differentially expressed
# 
# PARAMETERS:
#   phenoRow - single row of founders phenotype data
#   groupLabels - Specify which column of founders data belongs to group 0 and which to group 1.
#
# OUTPUT:
#  two p-values - for gene being up- and downregulated
#
############################################################################################################
findUsingTTest.internal <- function(phenoRow,groupLabels){
  a <- which(groupLabels==0)
  b <- which(groupLabels==1)
  if(mean(phenoRow[a],na.rm=T) < mean(phenoRow[b],na.rm=T)){
    what <- "less"
    return(c(0,t.test(phenoRow[a],phenoRow[b],alt=what)$p.value))
  }else{
    what <- "gre"
    return(c(t.test(phenoRow[a],phenoRow[b],alt=what)$p.value,0))
  }
}
