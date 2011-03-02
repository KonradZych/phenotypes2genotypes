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
	above <- 0
	below <- 0
	for(i in (data_vector)){
		if(i>=0.05){
			res <- res + 1
			above <- above + 1
		}else if(i<=-0.05){
			res <- res + 1
			below <- below +1
		}
	}
	output <- list(res,above,below)
	output
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
	colnames(geno_matrix)<-colnames(brassica_phenotypes, do.NULL = FALSE)
	rownames(geno_matrix)<-rownames(brassica_phenotypes, do.NULL = FALSE)
	geno_matrix
}

choose_right <- function(raw_expression_matrix,points_matrix,treshold,margin){
	output <- NULL
	for(i in 1:nrow(raw_expression_matrix)){
		non_zero <- points_matrix[i][[1]][[1]]
		above_min_below <- abs(points_matrix[i][[1]][[2]]-points_matrix[i][[1]][[3]])/2
		margin_range <- (non_zero*margin)
		if(non_zero<treshold){
			
		}else{
			if(above_min_below>margin_range){
			}else{
				output <- c(output,i)
			}
		}
	}
	result <- raw_expression_matrix[output,]
}

for( i in 1:nrow(qtl)){
for( j in 1:ncol(qtl)){
	if(qtl[i,j]==Inf){qtl[i,j]<-100}
}
}

setwd("D:/data")
library(basicQtl)
brassica <- as.matrix(read.table("Expression_BrassicaRapa_10chr2.txt",sep=""))
points_matrix <- select_diversed(brassica)
brassica_phenotypes <- choose_right(brassica, points_matrix,46,0.1)
brassica_genotypes <- create_genotypes(brassica_phenotypes)
brassica_genotypes_ord_1 <- un_neighbor(t(brassica_genotypes),1,500,10)
brassica_genotypes_ord_2 <- un_neighbor(brassica_genotypes_ord,1,5000,10)
brassica_genotypes_ord_3 <- un_neighbor(brassica_genotypes_ord,1,500,1)
brassica_genotypes_ord_4 <- un_neighbor(brassica_genotypes_ord_3,1,5000,10)
brassica_genotypes_ord <- un_neighbor(brassica_genotypes_ord,1,1000,10)
brassica_genotypes_cor <- cor((brassica_genotypes_ord),use="pairwise.complete.obs")
c_plot(selection[,10])