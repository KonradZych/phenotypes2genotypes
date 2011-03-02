chr_vect<-function(chr_vec){
res<-sum(apply(chr_vec,1,chr_vect2))
}

chr_vect2 <- function(char){
	if(substr(as.character(char),2,2)==substr(as.character(char),1,1)){
		return(2)
	}else{
		return(1)
	}
	}
	
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


diversity_in_row <- function(data_vector){
	res <- 0
	for(i in (data_vector)){
		if(i>=0.05||i<=-0.05){
			res <- res + 1
		}
	}
	res
}

select_diversed <- function(data_matrix){
	point_matrix <- apply(data_matrix,1,diversity_in_row)
	point_matrix
}

genotype_col <- function(brassica_phenotypes_vector){
	result <- as.vector(matrix(1,1,length(brassica_phenotypes_vector)))
	for(i in 1:length(brassica_phenotypes_vector)){
		if(brassica_phenotypes_vector[i]<0){
			result[i] <- 2
		}
	}
	result
}

create_genotypes <- function(brassica_phenotypes){
	geno_matrix <- apply(brassica_phenotypes,2,genotype_col)
	#geno_matrix <- matrix(geno_matrix,nrow(data_matrix),ncol(data_matrix))
	geno_matrix
}

for( i in 1:nrow(qtl)){
for( j in 1:ncol(qtl)){
	if(qtl[i,j]==Inf){qtl[i,j]<-100}
}
}

setwd("D:/data")
brassica <- as.matrix(read.table("Expression_BrassicaRapa_10chr2.txt",sep=""))
points_matrix <- select_diversed(brassica)
brassica_phenotypes <- brassica[which(points_matrix>=45),]
brassica_genotypes <- create_genotypes(brassica_phenotypes)
c_plot(selection[,10])