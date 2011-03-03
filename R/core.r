expressionScores_row <- function(expressionMatrix_row, no_expression_value=0, verbose=FALSE, debugMode=0){
	#Checkpoints
	if(!is.numeric(expressionMatrix_row)) stop("Data you provide to expressionScores_row in expressionMatrix_row must be numeric!")
	if(!is.numeric(no_expression_value)) stop("no_expression_value must be a numeric value")
	if(verbose!=TRUE && verbose!=FALSE) print("verbose must be boolean (TRUE/FALSE)")
	if(!pmatch(debugMode,c(0,1,2,3))) print("debugMode could be either 0,1,2 or 3")
	#Function itself
	if(debugMode==1) cat("ExpressionScores_row starting withour errors in checkpoint.")
	if(debugMode==2||debugMode==3) cat("ExpressionScores_row starting withour errors in checkpoint. Paramteres values - expressionMatrix_row:",expressionMatrix_row,"no_expression_value:",no_expression_value,"verbose:",verbose,"debugMode",debugMode,"\n")
	s <- proc.time()
	zeroes <- 0
	above <- 0
	below <- 0
	if(verbose) print("Starting loop.")
	for(i in 1:length(data_vector)){
		sl <- proc.time()
		if(i>no_expression_value){
			above <- above + 1
		}else if(i<no_expression_value){
			below <- below +1
		}else{
			res <- res + 1
		}
		el<-proc.time()
		if(verbose) cat("Done element:",i,"this took:",(el-sl)[3],"seconds,","left elements:",length(data_vector)-i,"estimated time to finish:",(length(data_vector)-i)*(el-s)[3]/i,"seconds\n.")
	}
	output <- list(res,above,below)
	e<-proc.time()
	if(verbose) cat("Done, expressionScores_row took:",(e-s)[3],"seconds\n.")
	if(debugMode==3) cat("expressionScores_row done, returns output:",output,"\n")
	output
}

#applies for every row function that counts 0, aboves and belows
expressionScores <- function(expressionMatrix, no_expression_value=0, verbose=FALSE, debugMode=0){
	#Checkpoints
	if(!is.numeric(expressionMatrix)) stop("Data you provide to expressionScores in expressionMatrix must be numeric!")
	if(!is.numeric(no_expression_value)) stop("no_expression_value must be a numeric value")
	if(verbose!=TRUE && verbose!=FALSE) print("verbose must be boolean (TRUE/FALSE)")
	if(!pmatch(debugMode,c(0,1,2,3))) print("debugMode could be either 0,1,2 or 3")
	if(no_expression_value<min(expressionMatrix)) stop("no_expression_value too small, lower than minimal value in matrix")
	#Function itself
	if(debugMode==1) cat("ExpressionScores_row starting withour errors in checkpoint.")
	if(debugMode==2||debugMode==3) cat("ExpressionScores_row starting withour errors in checkpoint. Paramteres values- expressionMatrix:",expressionMatrix,"no_expression_value",no_expression_value,"verbose:",verbose,"debugMode",debugMode,"\n")
	#Function itself
	s <- proc.time()
	output <- apply(expressionMatrix,1,expressionScores_row)
	e<-proc.time()
	if(verbose) cat("Done, expressionScores took:",(e-s)[3],"seconds\n.")
	if(debugMode==3) cat("expressionScores done, returns output:",output,"\n")
	output
}

#function that choses form the matrix only appropriate markers with specified rules
appriopiateMarkers <- function(expressionMatrix, proportion=50, margin=5, zeros_allowed = 0, no_expression_value=0, verbose=FALSE, debugMode=0){
	#Checkpoints
	if(!is.numeric(expressionMatrix)) stop("Data you provide to appriopiateMarkers in expressionMatrix must be numeric!")
	if(!is.numeric(proportion)) stop("proportion must be a numeric value")
	if(!is.numeric(margin)) stop("proportion must be a numeric value")
	if(!is.numeric(zeros_allowed)) stop("proportion must be a numeric value")
	if(!is.numeric(no_expression_value)) stop("no_expression_value must be a numeric value")
	if(verbose!=TRUE && verbose!=FALSE) print("verbose must be boolean (TRUE/FALSE)")
	if(!pmatch(debugMode,c(0,1,2,3))) print("debugMode could be either 0,1,2 or 3")	
	if(proportion < 1 || proportion > 99) stop("Proportion is a percentage (1,99)")
	if(zeros_allowed < 0 || zeros_allowed > ncol(expressionMatrix)) stop("Zeros_allowed is a number (0,lenght of the row)")
	if(margin < 0 || margin > proportion)) stop("Margin is a percentage (0,proportion)")
	if(no_expression_value<min(expressionMatrix)) stop("no_expression_value too small, lower than minimal value in matrix")
	if(debugMode==1) cat("appriopiateMarkers starting withour errors in checkpoint.")
	if(debugMode==2||debugMode==3) cat("appriopiateMarkers starting withour errors in checkpoint. Paramteres values- expressionMatrix:",expressionMatrix,"proportion",proportion,"margin",margin,"zeros_allowed",zeros_allowed,"no_expression_value",no_expression_value,"verbose:",verbose,"debugMode",debugMode,"\n")
	#Function itself
	s <- proc.time()
	#needs function providing points
	points_vector <- expressionScores(expressionMatrix)
	for(i in 1:nrow(expressionMatrix)){
		sl <- proc.time()
		non_zero <- points_vector[i][[1]][[1]]
		above_min_below <- abs(points_vector[i][[1]][[2]]-points_vector[i][[1]][[3]])/2
		margin_range <- (non_zero*margin)
		if(dim(expressionMatrix)[2]-non_zero < zeros_allowed){
			if(above_min_below>margin_range){
			}else{
				output <- c(output,i)
			}
		}
		el <- proc.time()
		if(verbose) cat("Done with element",i,"took:",(el-sl)[3],", estimated time remining:",(nrow(expressionMatrix)-i) * (el-sl)[3],"\n")
	}
	result <- expressionMatrix[output,]
	#we habe now vector which for every row in the matrix says us, how many zeros there are and how many above/below vaules
	e<-proc.time()
	if(verbose) cat("Done, expressionScores took:",(e-s)[3],"seconds\n.")
	if(debugMode==3) cat("expressionScores done, returns result:",result,"\n")
	result

}
