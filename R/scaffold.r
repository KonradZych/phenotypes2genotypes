#***************************************************************************#
#How should workflow look like:                                             #
# 1. raw data					                                            #
# 2. in c - handling big data, function to read column, row, dropping       #
#    markers that are evelny ex[ressed in most indyviduals					#





#utility
#c_plot - plotting routine, with three colors - three groups divided by cutoff
c_plot <- function(m_vector,cutoff=0,...){
	result <- as.vector(matrix("a",1,length(m_vector)))
	for(i in 1:length(m_vector)){
		if(m_vector[i]>cutoff){
		if(m_vector[i]>cutoff){
			result[i] <- "red"
		}else if(m_vector[i]<(-cutoff)){
			result[i] <- "blue"
		}else{
			result[i] <- "grey"
		}
	}
	plot(m_vector,col=result,...)
}
}



#internal
#diversity_in_row
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

#core
#select_diversed - uses diversity_in_row
select_diversed <- function(data_matrix,...){
	point_matrix <- apply(data_matrix,1,diversity_in_row,5,hello)
	point_matrix
}

#internal
#genotype_col
genotype_col <- function(brassica_phenotypes_vector){
	result <- as.vector(matrix(1,1,length(brassica_phenotypes_vector)))
	for(i in 1:length(brassica_phenotypes_vector)){
		if(brassica_phenotypes_vector[i]<0){
			result[i] <- 2
		}
	}
	result
}

#core
#create_genotypes - uses genotype_col
create_genotypes <- function(brassica_phenotypes){
	geno_matrix <- apply(brassica_phenotypes,2,genotype_col)
	#geno_matrix <- matrix(geno_matrix,nrow(data_matrix),ncol(data_matrix))
	colnames(geno_matrix)<-colnames(brassica_phenotypes, do.NULL = FALSE)
	rownames(geno_matrix)<-rownames(brassica_phenotypes, do.NULL = FALSE)
	geno_matrix
}

#core
#choose_right

choose_right(expression)

errormsg_proportion_tolow <- function(){
 return("Proportion is a percentage (1,99)")
}

choose_right <- function(expressionMatrix, proportion=50, margin=5, zeros_allowed = 0, verbose=FALSE, debugMode=0){
	if(proportion < 1) stop(errormsg_proportion_tolow())
	if(proportion > 99) stop("Proportion is a percentage (1,99)")
	output <- NULL
	s <- proc.time()
	for(i in 1:nrow(raw_expression_matrix)){
		sl <- proc.time()
		non_zero <- points_matrix[i][[1]][[1]]
		above_min_below <- abs(points_matrix[i][[1]][[2]]-points_matrix[i][[1]][[3]])/2
		margin_range <- (non_zero*margin)
		if(dim(raw_expression_matrix)[2]-non_zero < zeros_allowed)
			if(above_min_below>margin_range){
			}else{
				output <- c(output,i)
			}
		}
		e <- proc.time()
		if(verbose) cat("Done with element",i,"took:",(e-sl)[3],", estimated time remining:",(nrow(raw_expression_matrix)-i) * (e-sl)[3],"\n")
	}
	result <- raw_expression_matrix[output,]
	if(verbose) cat("Done with analysis:",(nrow(raw_expression_matrix)-i) * (e-s)[3]," Seconds\n")
	class(result) <- c(class(results),"expression_data")
}

setwd("D:/data")
library(basicQtl)
brassica <- as.matrix(read.table("Expression_BrassicaRapa_10chr2.txt",sep=""))
points_matrix <- select_diversed(brassica)
brassica_phenotypes <- choose_right(brassica, points_matrix,3,0.1)
brassica_genotypes <- create_genotypes(brassica_phenotypes)
brassica_genotypes_ord_1 <- un_neighbor(t(brassica_genotypes),1,500,10)
brassica_genotypes_ord_2 <- un_neighbor(brassica_genotypes_ord,1,5000,10)
brassica_genotypes_ord_2_cor <-cor(brassica_genotypes_ord_2,use="pairwise.complete.obs")
brassica_genotypes_ord_3 <- un_neighbor(brassica_genotypes_ord,1,500,1)
brassica_genotypes_ord_4 <- un_neighbor(brassica_genotypes_ord_3,1,5000,10)
brassica_genotypes_ord <- un_neighbor(brassica_genotypes_ord,1,1000,10)
brassica_genotypes_cor <- cor((brassica_genotypes_ord),use="pairwise.complete.obs")
c_plot(selection[,10])