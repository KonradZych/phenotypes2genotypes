#orderedCross - produces from genotypic matrix file containing object of type cross, reads it into R a returns
##Depends on read.cross (qtl library)
#chrom_matrix - matrix of genotypic data, rows - markers, cols - individuals
#nr_iterations - not yet used, maybe never will be
#groups - nr of groups we are dividing our dataset to, hopefully coresponds to nr of chromosomes
#outputFile - file where object of type cross is being saved
#verbose - standard
#debugMode - standard
#genos - argument passed to read.cross (chars describing genotypes)
#usage cross <- orderedCross(genotypicMatrix)
orderedCross <- function(chrom_matrix,nr_iterations=100,groups=10,outputFile="cross2.csv",verbose=FALSE,debugMode=0,genos=c("0","1")){
	library(qtl)
	#r <- un_best_clustering(chrom_matrix,nr_iterations,groups)
	if(verbose){cat("orderMarkers starting.\n\n")}
	s <- proc.time()
	
	#printing faked phenotype
	cat("",file=outputFile)
	cat("phenotype",sep="",file=outputFile,append=T)
	cat(",,,",sep="",file=outputFile,append=T)
	cat(paste(runif(48,0,100),collapse=","),sep="",file=outputFile,append=T)
	cat("\n",sep="",file=outputFile,append=T)
	
	#Dividing data in ten groups based on corelation using kmeans
	cor_matrix <- cor(t(chrom_matrix), use="pairwise.complete.obs")
	r <- kmeans(chrom_matrix,groups)
	sorted <- sort(r[[1]],index.return=T)

	#printing genotypic data to file in specified format
	for(i in 1:groups){
		sl <- proc.time()
		if(verbose){
			cat("   ### Writing chromosome: ",i,"nr of markers:",length(which(r[[1]]==i)),"###\n")
			cat(paste("      -> Marker ",colnames(sorted[[1]][which(sorted[[1]]==i)]),"\n",sep=""))
		}
		write.table(cbind(sorted[[1]][which(sorted[[1]]==i)], 1:table(sort(r[[1]]))[i], chrom_matrix[sorted[[2]][which(sorted[[1]]==i)],]) , file=outputFile , sep="," , quote=F , col.names=F , append=T)
		el <- proc.time()
		if(verbose){cat("   ### Writing chromosome:",i,"taken:",(el-sl)[3],"seconds. ###\n\n\n")}
	}

	
	#reading freshly made file to R
	cross <- read.cross("csvr",file=outputFile,genotypes=genos)
	#forcing cross time to RIL
	class(cross)[1] <- "riself"
	e <- proc.time()
	if(verbose){cat("Done without errors in:",(e-s)[3],"seconds.\n")}
	#returning cross
	invisible(cross)
}