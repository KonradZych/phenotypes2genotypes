crossParser <- function(genotypicMatrix,outputFile="cross.csv",verbose=FALSE, debugMode=0){
	genotypicMatrix <- switchMatrixValues(genotypicMatrix,"-",0,1,"A","B")
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

#switches values from NA, before1 and before2 to naValue, after1 and after2, respectively
switchMatrixValues <- function(matrix_to_be_cleaned,naValue=NA,before1,before2,after1,after2){
	matrix_to_be_cleaned[which(is.na(matrix_to_be_cleaned))]<-naValue
	matrix_to_be_cleaned[which(matrix_to_be_cleaned==before1)]<-after1
	matrix_to_be_cleaned[which(matrix_to_be_cleaned==before2)]<-after2
	matrix_to_be_cleaned
}


#recombinationCount - counting recobinations needen to go from one matrix to another
#genotypicMatrix - rows: markers, cols: individuals
recombinationCount <- function(genotypicMatrix){
	res <- apply(genotypicMatrix,1,recombinationCountRow,genotypicMatrix)
}