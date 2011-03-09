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

cEquals <- function(x,splitVal){
	sum(x==splitVal)
}

cLess <- function(x,splitVal){
	sum(x>splitVal)
}

cMore <- function(x,splitVal){
	sum(x<splitVal)
}

check <- function(x, overlapInd, margin_range ,splitVal){
	r <- FALSE
	if(cEquals(x,splitVal) <= overlapInd){
		if(abs(cMore(x,splitVal)-cLess(x,splitVal)) < margin_range){
			r <- TRUE
		}
	}
	r
}

#switches values from from1 and from2 to to1 and to2
switchMatrixValues <- function(Matrix,from1,from2,to1,to2){
	Matrix[which(Matrix==from1)] <- to1
	Matrix[which(Matrix==from2)] <- to2
	Matrix
}

f1 <- function(x,m){
	apply(m,1,f2,x)
}

f2 <- function(x,oldx){
	sum(x!=oldx)
}
 #recombination count, structure it!
res <- apply(brassica_genotypes,1,f1,brassica_genotypes)

#THIS COULD BE USEFUL
orderGroup <- function(chromMatrix,verbose=FALSE,debugMode=0){
	if(verbose)cat("   -> Starting un_order_chromosome_by_seriation\n",file="orderMarkers.log",append=T)
	cur <- abs(cor(chromMatrix,use="pairwise.complete.obs"))
	for(i in 1:10){
		si <- proc.time()
		if(verbose)cat("      -> Iteration",i,"\n",file="orderMarkers.log",append=T)
		o <- seriate(t((cur)))
		cur <- cur[get_order(o),get_order(o)]
		ei <- proc.time()
		if(verbose)cat("         -> Iteration",i,"taken:",(ei-si)[3],"seconds. Estimated time left:",(10-i)*(ei-si)[3],"seconds.\n",file="orderMarkers.log",append=T)
	}
	get_order(o)
}

brassica_genotypes2 <- matrix(as.numeric(switchMatrixValues(brassica_genotypes,"A","B",0,1)),nrow(brassica_genotypes),ncol(brassica_genotypes))


#chrom_matrix - matrix of genotypic data, rows - markers, cols - individuals
#nr_iterations
#groups
#outputFile
#verbose
#debugMode
#logFile
orderedCross <- function(chrom_matrix,nr_iterations=100,groups=10,outputFile="cross2.csv",verbose=FALSE,debugMode=0,logFile="orderMarkers.log"){
	
	#r <- un_best_clustering(chrom_matrix,nr_iterations,groups)
	cat("",file=logFile)
	if(verbose){cat("orderMarkers starting.\n\n",file=logFile,append=T)}
	s <- proc.time()
	cat("",file=outputFile)
	cat("phenotype",sep="",file=outputFile,append=T)
	cat(",,,",sep="",file=outputFile,append=T)
	cat(paste(runif(48,0,100),collapse=","),sep="",file=outputFile,append=T)
	cat("\n",sep="",file=outputFile,append=T)
	chrom_matrix <- matrix(as.numeric(switchMatrixValues(chrom_matrix,"A","B",0,1)),nrow(chrom_matrix),ncol(chrom_matrix))
	cor_matrix <- cor(t(chrom_matrix), use="pairwise.complete.obs")
	r <- kmeans(chrom_matrix,groups)
	sorted <- sort(r[[1]],index.return=T)
	#res <- NULL
	start <- 1
	for(i in 1:groups){
		sl <- proc.time()
		if(verbose){
			cat("   ### Writing chromosome: ",i,"nr of markers:",length(which(r[[1]]==i)),"###\n",file=logFile,append=T)
			cat(paste("      -> Marker ",colnames(sorted[[1]][which(sorted[[1]]==i)]),"\n",sep=""),file=logFile,append=T)
		}
		write.table(cbind(sorted[[1]][which(sorted[[1]]==i)], 1:table(sort(r[[1]]))[i], brassica_genotypes[sorted[[2]][which(sorted[[1]]==i)],]) , file=outputFile , sep="," , quote=F , col.names=F , append=T)
		el <- proc.time()
		if(verbose){cat("   ### Writing chromosome:",i,"taken:",(el-sl)[3],"seconds. ###\n\n\n",file=logFile,append=T)}
	}
	e <- proc.time()
	cross <- read.cross("csvr",file=outputFile,genotypes=c("A","B"))
	class(cross)[1] <- "riself"
	if(verbose){cat("orderMarkers taken:",(e-s)[3],"seconds.\n",file=logFile,append=T)}
	if(verbose){cat("Done without errors in:",(e-s)[3],"seconds. Saved log to:",logFile,"\n")}
	invisible(cross)
}



test.appriopriateMarkers <- function(){
	expressionMatrix <- as.matrix(read.table("Expression_BrassicaRapa_10chr2.txt",sep=""))
	brassica_genotypes <- appriopriateMarkers(expressionMatrix,margin=0.5, overlapInd=0, verb=T)
	dim(brassica_genotypes)
	cross <- orderedCross(brassica_genotypes,verbose=T)
	plot.rf(formLinkageGroups(cross,reorgMarkers=F))

}

#function that choses form the matrix only appropriate markers with specified rules
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
	if(proportion < 1 || proportion > 99) stop("Proportion is a percentage (1,99)")
	if(overlapInd < 0 || overlapInd > ncol(expressionMatrix)) stop("overlapInd is a number (0,lenght of the row).")
	if(margin < 0 || margin > proportion) stop("Margin is a percentage (0,proportion)")
	if(splitVal<min(expressionMatrix)) stop("splitVal too small, lower than minimal value in matrix")
	if(splitVal>max(expressionMatrix)) stop("splitVal too big, higher than maximal value in matrix")
	if(debugMode==1) cat("appriopriateMarkers starting withour errors in checkpoint.\n")
	if(debugMode==2||debugMode==3) cat("appriopriateMarkers starting withour errors in checkpoint. Paramteres values- expressionMatrix:",expressionMatrix,"proportion:",proportion,"margin",margin,"zeros_allowed:",zeros_allowed,"no_expression_value:",no_expression_value,"genotypeLabel1:",genotypeLabel1,"genotypeLabel2:",genotypeLabel2,"splitVal:",splitVal,"verbose:",verbose,"debugMode:",debugMode,"\n")

	#Selection of the probes matching to the specified parameters
	margin_range <- (ncol(expressionMatrix)*margin)/100
	suitedRows <- apply(expressionMatrix,1,check, overlapInd,margin_range,splitVal)#+splitValmargin)
	expressionMatrix <- expressionMatrix[suitedRows,]

	ep <- proc.time()
	if(verbose) cat("Selected proper probes, took:",(ep-s)[3],"seconds. Creating genotype matrix.\n")
	
	#Transform numeric values to genotypes
	expressionMatrix[which(expressionMatrix>splitVal)] <- genotypes[1]
	expressionMatrix[which(expressionMatrix<splitVal)] <- genotypes[2]
	expressionMatrix[which(expressionMatrix==splitVal)] <- NA

	eg <- proc.time()
	if(verbose) cat("Created genotype matrix, took:",(eg-ep)[3],"seconds.\n")
	e<-proc.time()
	if(verbose) cat("Done, appriopiateMarkers took:",(e-s)[3],"seconds.\n")

	expressionMatrix
}

#DON"T USE
cleanPlus<-function(matrix_to_be_cleaned,naValue,genoAValue,genoBValue){
	matrix_to_be_cleaned[which(is.na(matrix_to_be_cleaned))]<-naValue
	matrix_to_be_cleaned[which(matrix_to_be_cleaned==0)]<-genoAValue
	matrix_to_be_cleaned[which(matrix_to_be_cleaned==1)]<-genoBValue
	matrix_to_be_cleaned
}

crossParser <- function(genotypicMatrix,outputFile="cross.csv",verbose=FALSE, debugMode=0){
	genotypicMatrix <- cleanPlus(genotypicMatrix,"-","A","B")
	if(verbose){cat("crossParser: startin function.\n")}
	s <- proc.time()
	cat("",file=outputFile)
	cat("phenotype",sep="",file=outputFile,append=T)
	cat(",,,",sep="",file=outputFile,append=T)
	cat(paste(runif(48,0,100),collapse=","),sep="",file=outputFile,append=T)
	cat("\n",sep="",file=outputFile,append=T)
	e1 <- proc.time()
	if(verbose){cat("crossParser: started file, faked phenotype, strating genotypes. Taken:",(e1-s)[3],"s so far.\n")}
	if(verbose){cat("crossParser: starting loop, lines to write:",nrow(genotypicMatrix),"\n")}
	write.table(cbind(1,1:nrow(genotypicMatrix),genotypicMatrix),file=outputFile,sep=",",quote=F,col.names=F,append=T)
	e <- proc.time()
	if(verbose){cat("crossParser: finished. Taken:",(e-s)[3],"s. Result written to:",outputFile,"\n")}
	cross <- read.cross("csvr",file=outputFile,genotypes=c("A","B"))
	class(cross)[1] <- "riself"
	invisible(cross)
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
brassica_genotypes <- appriopriateMarkers("Expression_BrassicaRapa_10chr2.txt",zeros_allowed=0, no_expression_value=0.05)
brassica_genotypes_cor <- cor(t(brassica_genotypes), use="pairwise.complete.obs")
d <- dist(brassica_genotypes)
o <- seriate(d)
brassica_genotypes_cor <- cor(brassicaOrdered, use="pairwise.complete.obs")
brassica_genotypes_cor <- cor(t(brassica_genotypes)[,get_order(o)], use="pairwise.complete.obs")
appriopiateMarkers(m,zeros_allowed=0,verbose=TRUE,debugMode=3)
