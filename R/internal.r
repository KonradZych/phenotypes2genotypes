#INTERNAL
#WORK here
check_parameters <- function(functionName, numericParameters, verboseParameter, debugModeParameter, verbose=FALSE, debugMode=0){
	if(is.empty(verbose)) stop("Verbose must be boolean, but is empty rigth now!")
	if(is.empty(debugMode)) stop("DebugMode must be 0,1,2 or 3 but is empty rigth now!")
	if(verbose || debugMode==1 || debugMode==2 || debugMode==3) cat("Checking parameters for:",functionName,"\n")
	for(i in numericParameters){
		if(!is.numeric(i)||is.empty(i)) stop("Data you provide to ",functionName," must be numeric! And you provide it with: ", i)
	}
	if(verboseParameter!=TRUE && verboseParameter!=FALSE) stop("verbose must be boolean (TRUE/FALSE)")
	if(!pmatch(debugModeParameter,c(0,1,2,3))) stop("debugMode could be either 0,1,2 or 3")
}

#INTERNAL
cEquals <- function(x,splitVal){
	sum(x==splitVal)
}

#INTERNAL
cLess <- function(x,splitVal){
	sum(x>splitVal)
}

#INTERNAL
cMore <- function(x,splitVal){
	sum(x<splitVal)
}

#INTERNAL
check <- function(x, overlapInd, margin_range ,splitVal){
	r <- FALSE
	if(cEquals(x,splitVal) <= overlapInd){
		if(abs(cMore(x,splitVal)-cLess(x,splitVal)) < margin_range){
			r <- TRUE
		}
	}
	r
}

#INTERNAL
#recombinationCountRowSub - using recombinationCountRowSub to apply comparison between every two rows
recombinationCountRow <- function(genotypicMatrixRow,genotypicMatrix){
	apply(genotypicMatrix,1,recombinationCountRowSub,genotypicMatrixRow)
}

#INTERNAL
#recombinationCountRowSub - comparing two rows
recombinationCountRowSub <- function(genotypicMatrixRow1,genotypicMatrixRow2){
	sum(genotypicMatrixRow1!=genotypicMatrixRow2)
}