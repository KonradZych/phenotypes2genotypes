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

is.empty <- function(elementToBeChecked){
	if(is.na(elementToBeChecked&&1)||is.na(elementToBeChecked)){
	return(TRUE)
	}else{
	return(FALSE)
	}
}

expressionScores_row <- function(expressionMatrix_row, no_expression_value=0, verbose=FALSE, debugMode=0){
	#Checkpoints
	check_parameters("expressionScores_row",list(expressionMatrix_row,no_expression_value),verbose,debugMode,verbose,debugMode)
	#Function itself
	if(debugMode==1) cat("ExpressionScores_row starting withour errors in checkpoint.\n")
	if(debugMode==2||debugMode==3) cat("ExpressionScores_row starting withour errors in checkpoint. Paramteres values - expressionMatrix_row:",expressionMatrix_row,"no_expression_value:",no_expression_value,"verbose:",verbose,"debugMode",debugMode,"\n")
	s <- proc.time()
	zeroes <- 0
	above <- 0
	below <- 0
	if(debugMode==3) cat("Starting loop.\n")
	for(i in 1:length(expressionMatrix_row)){
		sl <- proc.time()
		if(expressionMatrix_row[i]>no_expression_value){
			above <- above + 1
		}else if(expressionMatrix_row[i]<(-no_expression_value)){
			below <- below + 1
		}else{
			zeroes <- zeroes + 1
		}
		el<-proc.time()
		if(debugMode==3) cat("Done element:",i,"this took:",(el-sl)[3],"seconds,","left elements:",length(expressionMatrix_row)-i,"estimated time to finish:",(length(expressionMatrix_row)-i)*(el-s)[3]/i,"seconds.\n")
	}
	output <- list(zeroes,above,below)
	e<-proc.time()
	if(verbose) cat("Done, expressionScores_row took:",(e-s)[3],"seconds.\n")
	if(debugMode==3){ 
		cat("expressionScores_row done, returns output:\n")
		print(output)
	}
	output
}

#applies for every row function that counts 0, aboves and belows
expressionScores <- function(expressionMatrix, no_expression_value=0, verbose=FALSE, debugMode=0){
	#Checkpoints
	check_parameters("expressionScores",list(expressionMatrix,no_expression_value),verbose,debugMode,verbose,debugMode)
	if(no_expression_value<min(expressionMatrix)) stop("no_expression_value too small, lower than minimal value in matrix")
	#Function itself
	if(debugMode==1) cat("expressionScores starting withour errors in checkpoint.\n")
	if(debugMode==2||debugMode==3) cat("expressionScores starting withour errors in checkpoint. Paramteres values- expressionMatrix:",expressionMatrix,"no_expression_value",no_expression_value,"verbose:",verbose,"debugMode",debugMode,"\n")
	#Function itself
	s <- proc.time()
	output <- apply(expressionMatrix,1,expressionScores_row,no_expression_value,verbose,debugMode)
	e<-proc.time()
	if(verbose) cat("Done, expressionScores took:",(e-s)[3],"seconds.\n")
	if(debugMode==3){ 
		cat("expressionScores done, returns output:\n")
		print(output)
	}
	output
}

#function that choses form the matrix only appropriate markers with specified rules
appriopriateMarkers <- function(fileName, proportion=50, margin=5, zeros_allowed = 3, no_expression_value=0, genotypeLabel1=0, genotypeLabel2=1, genotypeSplittingValue=0, verbose=FALSE, debugMode=0){
	#Checkpoints
	expressionMatrix <- as.matrix(read.table("Expression_BrassicaRapa_10chr2.txt",sep=""))
	check_parameters("appriopriateMarkers",list(expressionMatrix,no_expression_value,proportion,margin,zeros_allowed),verbose,debugMode,verbose,debugMode)
	if(proportion < 1 || proportion > 99) stop("Proportion is a percentage (1,99)")
	if(zeros_allowed < 0 || zeros_allowed > ncol(expressionMatrix)) stop("Zeros_allowed is a number (0,lenght of the row).")
	if(margin < 0 || margin > proportion) stop("Margin is a percentage (0,proportion)")
	if(no_expression_value<min(expressionMatrix)) stop("no_expression_value too small, lower than minimal value in matrix")
	if(no_expression_value>max(expressionMatrix)) stop("no_expression_value too big, higher than maximal value in matrix")
	if(debugMode==1) cat("appriopriateMarkers starting withour errors in checkpoint.\n")
	if(debugMode==2||debugMode==3) cat("appriopriateMarkers starting withour errors in checkpoint. Paramteres values- expressionMatrix:",expressionMatrix,"proportion:",proportion,"margin",margin,"zeros_allowed:",zeros_allowed,"no_expression_value:",no_expression_value,"genotypeLabel1:",genotypeLabel1,"genotypeLabel2:",genotypeLabel2,"genotypeSplittingValue:",genotypeSplittingValue,"verbose:",verbose,"debugMode:",debugMode,"\n")
	#Function itself
	s <- proc.time()
	#needs function providing points
	points_vector <- expressionScores(expressionMatrix,no_expression_value,verbose,debugMode)
	output <- NULL
	for(i in 1:nrow(expressionMatrix)){
		sl <- proc.time()
		zero <- points_vector[i][[1]][[1]]
		non_zero <- dim(expressionMatrix)[2]-zero
		above_min_below <- abs(points_vector[i][[1]][[2]]-points_vector[i][[1]][[3]])
		margin_range <- (non_zero*margin)/100
		if(zero <= zeros_allowed){
			if(above_min_below < margin_range){
				output <- c(output,i)
			}
		}
		el <- proc.time()
		if(verbose) cat("Done with element",i,"took:",(el-sl)[3],", estimated time remining:",(nrow(expressionMatrix)-i) * (el-sl)[3],"\n")
	}
	result <- expressionMatrix[output,]
	geno_matrix <- result
	ep <- proc.time()
	if(verbose) cat("Selected proper probes, took:",(ep-s)[3],"seconds. Creating genotype matrix.\n")
	geno_matrix[which(result<genotypeSplittingValue)] <- genotypeLabel1
	geno_matrix[which(result>genotypeSplittingValue)] <- genotypeLabel2
	geno_matrix[which(result==genotypeSplittingValue)] <- NA
	eg <- proc.time()
	if(verbose) cat("Created genotype matrix, took:",(eg-ep)[3],"seconds.\n")
	#we habe now vector which for every row in the matrix says us, how many zeros there are and how many above/below vaules
	e<-proc.time()
	if(verbose) cat("Done, appriopiateMarkers took:",(e-s)[3],"seconds.\n")
	if(is.empty(geno_matrix)) warning("Genotypic matrix for selected conditions is empty.\n")
	if(debugMode==3){ cat("appriopiateMarkers done, returns geno_matrix:\n")
	print(geno_matrix)}
	gc()
	geno_matrix	
}

cleanPlus<-function(matrix_to_be_cleaned){
	for(h in 1:nrow(matrix_to_be_cleaned)){
		for(w in 1:ncol(matrix_to_be_cleaned)){
			if(is.na(matrix_to_be_cleaned[h,w])){
				matrix_to_be_cleaned[h,w] <- as.character("")
			}else if(matrix_to_be_cleaned[h,w]==0){
				matrix_to_be_cleaned[h,w] <- as.character("A")
			}else{
				matrix_to_be_cleaned[h,w] <- as.character("B")
			}
		}
	}
	matrix_to_be_cleaned
}

crossParser <- function(genotypicMatrix){
	genotypicMatrix <- cleanPlus(genotypicMatrix)
	Tfile <- "output.csv"
	cat("",file=Tfile)
	cat("phenotype",sep="",file=Tfile,append=T)
	cat(",,,",sep="",file=Tfile,append=T)
	cat(paste(runif(48,0,100),collapse=","),sep="",file=Tfile,append=T)
	cat("\n",sep="",file=Tfile,append=T)
	for(i in 1:nrow(genotypicMatrix)){
		cat(rownames(genotypicMatrix)[i],1,i,sep=",",file = Tfile,append=T)
		for(w in 1:ncol(genotypicMatrix)-1){
			cat(genotypicMatrix[i,w],sep="",file = Tfile,append=T)
			cat(",",sep="",file=Tfile,append=T)
		}
		cat(genotypicMatrix[i,ncol(genotypicMatrix)],sep="",file = Tfile,append=T)
		cat("\n",file = Tfile,append=T)
	}
}


record_format_parser <- function(genotypicMatrix){
	genotypicMatrix <- cleanPlus(genotypicMatrix)
	Tfile <- "output.loc"
	cat("",file=Tfile)
	cat("name = super\n",file = Tfile,append=T)
	cat("popt = RI8\n",file = Tfile,append=T)
	cat("nloc =",ncol(genotypicMatrix),"\n",file = Tfile,append=T)
	cat("nind =",nrow(genotypicMatrix),"\n\n",file = Tfile,append=T)
	for(i in 1:nrow(genotypicMatrix)){
		cat(rownames(genotypicMatrix)[i],"\n",file = Tfile,append=T)
		len <- ncol(genotypicMatrix)
		while(len>=50){
			cat("  ",file = Tfile,append=T)
			for(j in 1:10){
				cat(genotypicMatrix[i,(len-4):len],sep="",file = Tfile,append=T)
				cat(" ",file = Tfile,append=T)
				len <- len-5
			}
			cat("\n",file = Tfile,append=T)
		}
		cat("  ",file = Tfile, append=T)
		while(len>=5){
			cat(genotypicMatrix[i,(len-4):len],sep="",file = Tfile,append=T)
			cat(" ",file = Tfile,append=T)
			len <- len-5
		}
		for(w in 1:len){
			cat(genotypicMatrix[i,w],sep="",file = Tfile,append=T)
		}
		cat("\n",file = Tfile,append=T)
	}
}
record_format_parser(brassica_genotypes)

setwd("D:/data")
library(basicQtl)
brassica_genotypes <- appriopriateMarkers("Expression_BrassicaRapa_10chr2.txt",zeros_allowed=3, no_expression_value=0.05)
brassica_genotypes_cor <- cor(t(brassica_genotypes), use="pairwise.complete.obs")
d <- dist(brassica_genotypes)
o <- seriate(d)
brassica_genotypes_cor <- cor((un_ord), use="pairwise.complete.obs")
brassica_genotypes_cor <- cor(t(brassica_genotypes)[,get_order(o)], use="pairwise.complete.obs")
appriopiateMarkers(m,zeros_allowed=0,verbose=TRUE,debugMode=3)
