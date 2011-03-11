#WORK here
check_parameters <- function(functionName, numericParameters, booleanParameters, debugModeParameters, verbose=FALSE, debugMode=0){
	if(verbose || debugMode==1 || debugMode==2 || debugMode==3) cat("Checking parameters for:",functionName,"\n")
	if(numericParameters){if(!lapply(numericParameters,is.numeric))stop("One of parameters you provided to",functionName,"is not numeric. Check help file.\n")}
	if(booleanParameters){if(!lapply(booleanParameters,))stop("Verbose you provide to",functionName,"must be boolean. Check help file.\n")}
	if(debugModeParameters){if(!lapply(debugModeParameters,))stop("DebugMode you provide to",functionName,"must be 0,1,2 or 3. Check help file.\n")}
}

isBoolean <- function(x){
	if(x!=TRUE&&x!=FALSE){ return(FALSE)}
	else{ return(TRUE) }
}

isDebug <- function(x){
	if(!pmatch(x,c(0,1,2,3))){ return(FALSE)}
	else{ return(TRUE) }
}

cEquals <- function(x,splitVal){
	sum(x==splitVal)
}

cLess <- function(x,splitVal){
	sum(x>splitVal)
}

cMore <- function(x,splitVal){
	sum(x<splitVal)
}

#x
#overlapInd
#
check <- function(x, overlapInd, margin_range ,splitVal){
	r <- FALSE
	if(cEquals(x,splitVal) <= overlapInd){
		if(abs(cMore(x,splitVal)-cLess(x,splitVal)) < margin_range){
			r <- TRUE
		}
	}
	r
}

#recombinationCountRowSub - using recombinationCountRowSub to apply comparison between every two rows
recombinationCountRow <- function(genotypicMatrixRow,genotypicMatrix,flip=0,verbose=FALSE,debugMode=0){
	if(flip==0){output <- apply(genotypicMatrix,1,recombinationCountRowSub,genotypicMatrixRow)}
	if(flip==1){output <- apply(genotypicMatrix,1,recombinationCountRowFlipSub,genotypicMatrixRow)}
	output
}

#recombinationCountRowSub - comparing two rows
recombinationCountRowSub <- function(genotypicMatrixRow1,genotypicMatrixRow2){
	sum((genotypicMatrixRow1)!=(genotypicMatrixRow2))
}

#recombinationCountRowSub - comparing two rows
recombinationCountRowFlipSub <- function(genotypicMatrixRow1,genotypicMatrixRow2){
	sum((genotypicMatrixRow1)!=(1-(genotypicMatrixRow2)))
}

flipMarkers <- function(genotypicMatrix){
	brassicaReco <- recombinationCount(genotypicMatrix)
	brassicaRecoflipped <- recombinationCount(genotypicMatrix,flip=1)
	brassicaRecorows <- apply(brassicaReco,1,mean)
	brassicaRecoflippedrows <- apply(brassicaRecoflipped,1,mean)
	result <- NULL
	for(i in 1:length(brassicaRecoflippedrows)){
		cat(i,"\n")
		if(brassicaRecorows[i]<brassicaRecoflippedrows[i]){
			result <- rbind(result, brassicaReco[i,])
		}else{
			result <- rbind(result, brassicaRecoflipped[i,])
		}
	}
	result
}