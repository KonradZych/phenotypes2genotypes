#c_plot - no documantation yet!!!
c_plot <- function(m_vector,cutoff=0){
	result <- as.vector(matrix("a",1,length(m_vector)))
	for(i in 1:length(m_vector)){
		if(m_vector[i]>cutoff){
			result[i] <- "red"
		}else if(m_vector[i]<(-cutoff)){
			result[i] <- "blue"
		}else{
			result[i] <- "grey"
		}
	}
	plot(m_vector,col=result)
}

flip_plot <- function(vec,vecf){
	result <- as.vector(matrix("white",1,length(vec)))
	for(i in 1:length(vec)){
		if(vec[i]>vecf[i]){
			result[i] <- "red"
		}
	}
	plot(vecf,col=result)
}

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
switchMatrixValues <- function(matrix_to_be_cleaned,naValue=NA,before,after){
	matrix_to_be_cleaned[which(is.na(matrix_to_be_cleaned))]<-naValue
	matrix_to_be_cleaned[which(matrix_to_be_cleaned==before[1])]<-after[1]
	matrix_to_be_cleaned[which(matrix_to_be_cleaned==before[2])]<-after[2]
	matrix_to_be_cleaned
}


#recombinationCount - counting recobinations needen to go from one matrix to another
#genotypicMatrix - rows: markers, cols: individuals
recombinationCount <- function(genotypicMatrix,flip=0,verbose=FALSE,debugMode=0){
	genotypicMatrix <- switchMatrixValues(genotypicMatrix,before=c("A","B"),after=c(0,1))
	res <- apply(genotypicMatrix,1,recombinationCountRow,genotypicMatrix,flip)
	res
}

#c_plot - no documantation yet!!!
c_plot <- function(m_vector,cutoff=0){
	result <- as.vector(matrix("a",1,length(m_vector)))
	for(i in 1:length(m_vector)){
		if(m_vector[i]>cutoff){
			result[i] <- "red"
		}else if(m_vector[i]<(cutoff)){
			result[i] <- "blue"
		}else{
			result[i] <- "grey"
		}
	}
	plot(m_vector,col=result)
}
