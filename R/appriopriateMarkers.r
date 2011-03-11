#function that choses form the matrix only appropriate markers with specified rules
#Depends on: checkm (internal.R) 
#expressionMatrix -> columns -. indivudlas,rows markers
#proportion -> Proportion of individuals expected to carrying a certain genotype
#margin -> Proportion is allowed to varry between this margin (2 sided)
#splitVal -> Mean expression at which we split the indivuals by
#overlapInd -> number of individuals that are allowed to overlap (have splitVal)
#genotypes -> User defined genotypes for the output matrix
#verbose standrad
#debugmode standrad
appriopriateMarkers <- function(expressionMatrix, proportion = 50, margin = 5, splitVal = 0, splitValmargin = 0, overlapInd = 0, genotypes = c("A","B"), verbose=FALSE, debugMode=0){
	s <- proc.time()
	#source("internal.r") #We need check
	if(proportion < 1 || proportion > 99) stop("Proportion is a percentage (1,99)")
	if(overlapInd < 0 || overlapInd > ncol(expressionMatrix)) stop("overlapInd is a number (0,lenght of the row).")
	if(margin < 0 || margin > proportion) stop("Margin is a percentage (0,proportion)")
	if(splitVal<min(expressionMatrix)) stop("splitVal too small, lower than minimal value in matrix")
	if(splitVal>max(expressionMatrix)) stop("splitVal too big, higher than maximal value in matrix")
	if(debugMode==1) cat("appriopriateMarkers starting withour errors in checkpoint.\n")
	if(debugMode==2||debugMode==3) cat("appriopriateMarkers starting withour errors in checkpoint. Paramteres values- expressionMatrix:",expressionMatrix,"proportion:",proportion,"margin",margin,"zeros_allowed:",zeros_allowed,"no_expression_value:",no_expression_value,"genotypeLabel1:",genotypeLabel1,"genotypeLabel2:",genotypeLabel2,"splitVal:",splitVal,"verbose:",verbose,"debugMode:",debugMode,"\n")

	#Selection of the probes matching to the specified parameters
	margin_range <- (ncol(expressionMatrix)*margin)/100
	suitedRows <- apply(expressionMatrix,1,check,overlapInd,margin_range,splitVal)#+splitValmargin)
	expressionMatrix <- expressionMatrix[suitedRows,]

	ep <- proc.time()
	if(verbose) cat("Selected proper probes, took:",(ep-s)[3],"seconds. Creating genotype matrix.\n")
	
	#Transform numeric values to genotypes
	expressionMatrix2 <- expressionMatrix
	expressionMatrix2[which(expressionMatrix>splitVal)] <- genotypes[1]
	expressionMatrix2[which(expressionMatrix<splitVal)] <- genotypes[2]
	expressionMatrix2[which(expressionMatrix==splitVal)] <- NA

	eg <- proc.time()
	if(verbose) cat("Created genotype matrix, took:",(eg-ep)[3],"seconds.\n")
	e<-proc.time()
	if(verbose) cat("Done, appriopiateMarkers took:",(e-s)[3],"seconds.\n")

	expressionMatrix2
}

#counts how many times x eqauls splitVal
cEquals <- function(x,splitVal){
	sum(x==splitVal)
}

#counts how many times x is less than splitVal
cLess <- function(x,splitVal){
	sum(x>splitVal)
}

#counts how many times x is more than splitVal
cMore <- function(x,splitVal){
	sum(x<splitVal)
}

#check - checks if x meets specified requirments:
#splitVal - numeric value used to split the data into three subsets - eqauls, less and more
#overlapInd - how many individuals could be overalpping -> how many elememnts of x are eqaul to splitVal
#margin_range - we assume that nr of less and more should be the same, at least difference should be less than margin_range
check <- function(x, splitVal, overlapInd, margin_range){
	r <- FALSE
	if(cEquals(x,splitVal) <= overlapInd){
		if(abs(cMore(x,splitVal)-cLess(x,splitVal)) < margin_range){
			r <- TRUE
		}
	}
	r
}

#copy-paste into R for analysis
workflow.appriopriateMarkers <- function(){
	setwd("D:/data")
	library(basicQtl)
	library(qtl)
	expressionMatrix <- as.matrix(read.table("Expression_BrassicaRapa_10chr2.txt",sep=""))
	brassica_genotypes <- appriopriateMarkers(expressionMatrix,margin=0.5,genotypes=c(1,0),overlapInd=0, verb=T)
	brassicaFlipped <- flipMarkers(brassica_genotypes)
	cross <- orderedCross(brassica_genotypes,verbose=T)
	plot.rf(formLinkageGroups(cross,reorgMarkers=F))
	brassica_genotypes <- appriopriateMarkers(expressionMatrix,margin=0.5, overlapInd=0, verb=T)
	brassica_genotypes <- switchMatrixValues(brassica_genotypes,before=c("A","B"),after=c(0,1))
	brassicaReco <- recombinationCount(brassica_genotypes)
	brassicaRecoflipped <- recombinationCount(brassica_genotypes,flip=1)
}

#test function, doesn't need documenation
test.appriopriateMarkers <- function(verbose=FALSE, debugMode=0){
	correctData <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1),1,48)
	setwd("D:/data")
	library(basicQtl)
	library(qtl)
	expressionMatrix <- as.matrix(read.table("Expression_BrassicaRapa_10chr2.txt",sep=""))
	brassica_genotypes <- appriopriateMarkers(expressionMatrix, margin=0.5,genotypes=c(1,0),overlapInd=0, verbose=verbose, debugMode=debugMode)
	if(sum(which(brassica_genotypes[1,]!=correctData))){
		if(debugMode!=3){
			cat("Test went wrong! Starting test in debugMode=3, verbose=TRUE\n")
			test.appriopriateMarkers(verbose=TRUE, debugMode=3)
		}else{
			stop("Done test with debugMode=3, went wrong, stopping!")
		}
	}else{
		cat("Everything ok, appriopriateMarkers ready to serve, master!\n")
	}	
}