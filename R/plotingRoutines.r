##################################################################################################
#
# plottingRoutines.R
#
# Copyright (c) 2011, Konrad Zych
#
# Modified by Danny Arends
# 
# first written March 2011
# last modified July 2011
# last modified in version: 0.8.5
# in current version: active, not in main workflow
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
#
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
#
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: plotParentalExpression, plotChildrenExpression, plotMapComparison, plotMarkerDistribution
#				getChromosome.internal, getYLocs.internal, makeChromPal.internal,
#				makeTransPal.internal, switchChromosomes.internal
#           
#
##################################################################################################

############################################################################################################
#									*** plotParentalExpression ***
#
# DESCRIPTION:
#	plot red points for expression values for parent of type 0, blue for parent 1 and green lines for means
#	of rows
# 
# PARAMETERS:
# 	population - Ril type object, must contain founders phenotypic data.
# 	markers - markers to be printed numbers or names 
# 	groupLabels - Specify which column of founders data belongs to group 0 and which to group 1
# 
# OUTPUT:
#	plot
#
############################################################################################################
plotParentalExpression <- function(population, markers=1:100, groupLabels=c(0,0,1,1)){
	#*******checks*******
	is.population(population)
	
	#*******remove too short chromosomes*******
	if(nrow(population$founders$phenotypes)<length(markers)){
		cat("WARNING: There are only",nrow(population$founders$phenotypes),"markers in population$founders$phenotypes, function will plot them all.\n")
		markers <- 1:nrow(population$founders$phenotypes)
	}
	for(i in markers) if(!(i %in% 1:nrow(population$founders$phenotypes))){
		stop("ERROR: There is no marker number: ",i,", stopping \n")
	}
	founders <- population$founders$phenotypes[markers,]
	plot(x=markers[1], y=founders[1,1], xlim=c(min(markers),max(markers)), ylim=c(min(founders),max(founders)), col="red",
	xlab="Marker", ylab="Expression value", main="Parental gene expression data")
	for(i in 1:nrow(founders)){
		for(j in which(groupLabels==0)){
			if(i%%50==0) cat(i,j,founders[i,j],"\n")
			points(x=markers[i],y=founders[i,j],col="red")
		}
		for(k in which(groupLabels==1)){
			if(i%%50==0) cat(i,k,founders[i,k],"\n")
			points(x=markers[i],y=founders[i,k],col="blue")
		}
	}
	points(x=markers,y=apply(founders,1,mean),col="green", pch=95, cex=3)
}

############################################################################################################
#									*** plotChildrenExpression ***
#
# DESCRIPTION:
#	boxplot of data for selected markers + points of founders mean for each marker
# 
# PARAMETERS:
# 	population - Ril type object, must contain founders phenotypic data.
# 	markers - markers to be printed numbers or names 
# 
# OUTPUT:
#	plot
#
############################################################################################################
plotChildrenExpression <- function(population, markers=1:100){
	### checks
	is.population(population)
	
	### function itself
	if(nrow(population$offspring$phenotypes)<length(markers)){
		cat("WARNING: There are only",nrow(population$offspring$phenotypes),"markers in population$offspring$phenotypes, function will plot them all.\n")
		markers <- 1:nrow(population$offspring$phenotypes)
	}
	if(nrow(population$founders$phenotypes)<length(markers)){
		cat("WARNING: There are only",nrow(population$founders$phenotypes),"markers in population$founders$phenotypes, function will plot them all.\n")
		markers <- 1:nrow(population$founders$phenotypes)
	}
	for(i in markers) if(!(i %in% 1:nrow(population$offspring$phenotypes))){
		stop("ERROR: There is no marker number: ",i,", stopping \n")
	}
	for(i in markers) if(!(i %in% 1:nrow(population$founders$phenotypes))){
		stop("ERROR: There is no marker number: ",i,", stopping \n")
	}
	children <- population$offspring$phenotypes[markers,]
	founders <- population$founders$phenotypes[which(rownames(population$founders$phenotypes) %in% rownames(children)),]
	rownames(children) <- markers
	boxplot(t(children), ylim=c(min(children), max(children)),	xlab="Marker", ylab="Expression value", main="Children gene expression data")
	points(apply(founders,1,mean),col="green", pch=95, cex=3)
	points(apply(founders,1,max),col="red", pch=24, cex=1)
	points(apply(founders,1,min),col="blue", pch=25, cex=1)
}

############################################################################################################
#									*** plotMapComparison ***
#
# DESCRIPTION:
#	boxplot of data for selected markers + points of founders mean for each marker
# 
# PARAMETERS:
# 	cross - object of R/qtl cross type
# 	coloringMode - 1 - rainbow colors 2 - black for cis and red for trans located markers
# 	map - which map should be used for comparison:
#			- genetic - genetic map from cross$maps$genetic
#			- physical - physical map from cross$maps$physical
# 
# OUTPUT:
#	plot
#
############################################################################################################
plotMapComparison <- function(cross, coloringMode=1,map=c("genetic","physical")){ 
	
	crossContainsMap.internal(cross, map)
	#*******objects containing all information needen for function execution*******
	ys <- getYLocs.internal(cross)
	if(map=="genetic"){
			xs <- cross$maps$genetic[rownames(ys[[1]]),]
			#*******chromosomes lengths*******
			referenceChrom <- chromosomesLengths.internal(cross$maps$genetic)
			xs[,2] <- xs[,2] + referenceChrom[xs[,1]]
	}else if(map=="physical"){
			xs <- cross$maps$physical[rownames(ys[[1]]),]
			#*******chromosomes lengths*******
			referenceChrom <- chromosomesLengths.internal(cross$maps$physical)
			xs[,2] <- xs[,2] + referenceChrom[xs[,1]]
	}
	
	#*******positions of markers*******
	predictedLocs <- ys[[1]][,-1]
	referenceLocs <- xs[,2]
	
	#*******chromosomes lengths*******
	predictedChrom <- ys[[2]]
	
	#*******chromosome labels*******
	predictedChromLabels <- names(table(ys[[1]][,1]))
	referenceChromLabels <- names(table(xs[,1]))
	
	#*******chromosome labels positions*******
	predictedChromPos <- vector(mode="numeric",length(predictedChrom)-1)
	for(i in 1:length(predictedChrom)-1){
		predictedChromPos[i] <- (predictedChrom[i] + predictedChrom[i+1])/2
	}
	predictedChromPos[length(predictedChrom)] <- (predictedChrom[length(predictedChrom)] + max(predictedLocs))/2
	
	referenceChromPos <- vector(mode="numeric",length(referenceChrom)-1)
	for(i in 1:length(referenceChrom)-1){
		referenceChromPos[i] <- (referenceChrom[i] + referenceChrom[i+1])/2
	}
	#referenceChromPos[length(referenceChrom)] <- (referenceChrom[length(referenceChrom)] + max(referenceLocs))/2
	#*******color palette*******
	if(coloringMode==1){ 
		color <- makeChromPal.internal(ys[[1]],xs)
	}else if(coloringMode==2){
		color <- makeTransPal.internal(ys[[1]],xs)
	}
	
	#*******results of lin regr for each chromosome*******
	l <- vector(mode="list",length(table(ys[[1]][,1])))
	for(i in 1:length(table(ys[[1]][,1]))){
		a <- ys[[1]][which(ys[[1]][,1]==i),-1]
		b <- xs[which(rownames(xs) %in% names(a)),2]
		l[[i]] <- lm(a~b)$coefficients 
		p <- anova(lm(a~b))[[5]][1]
		cat("Chromosome",i,"lr coefficients",l[[i]],"corected by length diff",l[[i]][2]*((max(b)-min(b))/(max(a)-min(a))),"p-val",p,"\n")
	}
	#*******plotting points*******
	plot(x=referenceLocs, y=predictedLocs, xlim=c(min(referenceLocs),max(referenceLocs)), ylim=c(min(predictedLocs),max(predictedLocs)),
		xaxt="n", yaxt="n", col=color[[1]], pch=color[[2]], xlab="Reference map", ylab="Predicted map", main="Comparison of genetic maps")
	
	#*******adding chromosome labels and tics*******
	axis(1, at = referenceChrom[-1],labels = FALSE)
	axis(1, at = referenceChromPos,labels = referenceChromLabels, lwd = 0, tick = FALSE)
	axis(2, at = predictedChrom[-1],labels = FALSE)
	axis(2, at = predictedChromPos,labels = predictedChromLabels, lwd = 0, tick = FALSE)
	
	#*******adding marker tics*******
	axis(1, at = referenceLocs,labels = FALSE)
	#if(!is.null(cross$removed)) axis(1, at = cross$removed[,2],labels = FALSE, col.ticks = "red")
	axis(2, at = predictedLocs,labels = FALSE)
	
	#*******adding lines marking chromosomes ends*******
	for(x in 2:length(referenceChrom)){
		abline(v=sum(referenceChrom[x]),lty=2)
	}
	for(y in 2:length(predictedChrom)){
		abline(h=predictedChrom[y],lty=2)
	}
}

############################################################################################################
#									*** getChromosome.internal ***
#
# DESCRIPTION:
#	subfunction of plotMapComparison, returning list of chromosome numbers for all markers in cross
# 
# PARAMETERS:
# 	cross - object of class cross 
# 
# OUTPUT:
#	vector of numbers
#
############################################################################################################
getChromosome.internal <- function(cross){
	invisible(rep(1:length(nmar(cross)),nmar(cross)))
}

############################################################################################################
#									*** getYLocs.internal***
#
# DESCRIPTION:
#	subfunction of plotMapComparison, returning list of location of all markers in cross
# 
# PARAMETERS:
# 	cross - object of class cross 
#
# OUTPUT:
#	vector of numbers
#
############################################################################################################
getYLocs.internal <- function(cross){
	locs <- lapply(cross$geno,function(x){as.numeric(x$map)})
	chrlength <- as.vector(unlist(lapply(locs,max)),mode="numeric")
	locs <- as.numeric(unlist(locs))
	summaryLengths <- rep(0,length(chrlength))
	for(x in 2:length(chrlength)){
		summaryLengths[x] <- max(chrlength[x-1]) + summaryLengths[x-1] + 0.15 * max((chrlength))
	}
	chrids <- getChromosome.internal(cross)
	result <- matrix(0,length(locs),2)
	result[,1] <- chrids
	result[,2] <- summaryLengths[chrids]+locs
	rownames(result) <- markernames(cross)
	invisible(list(result,summaryLengths))
}

############################################################################################################
#									*** makeChromPal.internal ***
#
# DESCRIPTION:
#	subfunction of plotMapComparison, returning color pallete (rainbow colors)
# 
# PARAMETERS:
# 	ys1 - object of plotMapComparison, containing info about predicted map
# 	xs - object of plotMapComparison, containing info about reference map
#
# OUTPUT:
#	list of vector of colors (characters) and vector of numbers (symbol identifiers)
#
############################################################################################################
makeChromPal.internal <- function(ys1,xs){
	color <- vector(mode="character",nrow(ys1))
	names(color) <- rownames(ys1)
	symbol <- vector(mode="numeric",nrow(xs))
	names(symbol) <- rownames(xs)
	cl <- rainbow(length(table(ys1[,1])))
	for(i in rownames(ys1)){
		color[i] <- cl[ys1[i,1]]
		symbol[i] <- xs[i,1]
	}
	invisible(list(color, symbol))
}

############################################################################################################
#									*** makeTransPal.internal ***
#
# DESCRIPTION:
#	subfunction of plotMapComparison,returning color pallete (red/black)
# 
# PARAMETERS:
# 	ys1 - object of plotMapComparison, containing info about predicted map
# 	xs - object of plotMapComparison, containing info about reference map
#
# OUTPUT:
#	list of vector of colors (characters) and vector of numbers (symbol identifiers)
#
############################################################################################################
makeTransPal.internal <- function(ys1,xs){
	color <- vector(mode="character",nrow(ys1))
	names(color) <- rownames(ys1)
	symbol <- vector(mode="numeric",nrow(xs))
	names(symbol) <- rownames(xs)
	for(i in rownames(ys1)){
		if(ys1[i,1]==xs[i,1]){
			color[i] <- "black"
		}else{
			color[i] <- "red"
		}
		
		symbol[i] <- 1
	}
	invisible(list(color, symbol))
}

############################################################################################################
#									*** chromosomesLengths.internal ***
#
# DESCRIPTION:
#	function calculating lengths of chromosomes in given map
# 
# PARAMETERS:
# 	map - genetic or physical map (matrix with two cols - 1 - chromome nr, 2- postion on chromosome),
#		rownames - names of markers
#
# OUTPUT:
#	vector of lengths of chromosomes
#
############################################################################################################
chromosomesLengths.internal <- function(map){
	lengths <- vector(mode="numeric",length=(max(map[,1])+1))
	lengths[1] <- 0
	for(i in 1:max(map[,1])){
		lengths[i+1] <- max(map[which(map[,1]==i),2]) + lengths[i]
	}
	invisible(lengths)
}

############################################################################################################
#									*** plotMarkerDistribution ***
#
# DESCRIPTION:
#	plotting histogram of distribution of values for single marker and specified number of normal distribution 
#	curves, fitted to data using EM algorithm
# 
# PARAMETERS:
# 	population - population type object, must contain founders phenotypic data.
# 	marker - marker to be printed numbers or names 
# 	nrDistributions - numbers of normal distributions to be fitted
# 	logarithmic - TRUE - log(data) will be used instead of raw data
#
# OUTPUT:
#	plot
#
############################################################################################################
plotMarkerDistribution <- function(population,marker,nrDistributions,logarithmic=FALSE){
	if(missing(population)) stop("population object not provided.\n")
	if(missing(marker)) stop("marker not specified.\n")
	if(missing(nrDistributions)) stop("nrDistributions not specified.\n")
	if(length(marker)!=1) stop("plotMarkerDistribution can plot only one marker art once.\n")
	is.population(population)
	phenotypeRow <- population$offspring$phenotypes[marker,]
	if(logarithmic) phenotypeRow <- log(phenotypeRow)
	EM<-normalmixEM(phenotypeRow, k=nrDistributions)
	if(logarithmic){
		xlab <- "log(expression values)"
	}else{
		xlab <- "expression values"
	}
	h <- NULL
	len <- vector(mode="numeric",length=nrDistributions)
	for(i in 1:nrDistributions){
		len[i]<-length(phenotypeRow)*EM$lambda[which(EM$mu==sort(EM$mu)[i])]
		startVal <- sum(len[1:i-1])
		x <- sort(phenotypeRow)[startVal:(startVal+len[i])]
		h[i] <- hist(x,breaks=50)
	}
 	h0<-hist(phenotypeRow,breaks=50,col="grey62",border="grey70",xlab=xlab,ylab="Number of counts",main="Distribution of expression values for selected marker")
	colorP <- vector(mode="character",length=3)
	colorP[1] <- rgb(1,0,0)
	colorP[2] <- rgb(0,1,0)
	colorP[3] <- rgb(0,0,1)
	for(i in 1:nrDistributions){
		abline(v=EM$mu[i],col=colorP[i%%3+1])
		abline(v=c(EM$mu[i]-EM$sigma[i],EM$mu[i]+EM$sigma[i]),col=colorP[i%%3+1],lty=2)
		startVal <- sum(len[1:i-1])
		x <- sort(phenotypeRow)[startVal:(startVal+len[i])]
		xfit<-seq(min(x),max(x),length=40) 
		yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
		yfit <- yfit*diff(h[[i]][1:2])*length(x) 
		lines(xfit, yfit, col=colorP[i%%3+1], lwd=2)
	}
}

############################################################################################################
#									*** chromosomesCorelationPlot ***
#
# DESCRIPTION:
#	plot image show correlation between chromosomes of maps inside cross object
# 
# PARAMETERS:
# 	cross - object of R/qtl cross type
# 	map - which map should be used for comparison:
#			- genetic - genetic map from cross$maps$genetic
#			- physical - physical map from cross$maps$physical
# 
# OUTPUT:
#	plot
#
############################################################################################################
chromosomesCorelationPlot <- function(cross, map=c("genetic","physical")){ 
	inListCheck.internal(map,"map",c("genetic","physical"))
	crossContainsMap.internal(cross,map)
	if(map=="genetic"){
		cur_map <- cross$maps$genetic
	}else if(map=="physical"){
		cur_map <- cross$maps$physical
	}
	image(chromosomesCorelationMatrix.internal(cross,cur_map))
}

############################################################################################################
#									*** chromosomesCorelationMatrix.internal ***
#
# DESCRIPTION:
#	function calculating matrix of corelations between chromosomes
# 
# PARAMETERS:
# 	cross - object of R/qtl cross type
# 	cur_map - map to be used in comparison
# 
# OUTPUT:
#	matrix of corelations
#
############################################################################################################
chromosomesCorelationMatrix.internal <- function(cross, cur_map){
	### add mapMarkers here and in bestCorelated
	knchrom <- length(table(cur_map[,1]))
	result <- matrix(0, sum(nmar(cross)), nrow(cur_map))
	rownames(result) <- markernames(cross)
	colnames(result) <- rownames(cur_map)
	for(i in 1:length(cross$geno)){
		cur_ys <- cross$geno[[i]]$data[,]
		for(j in 1:knchrom){
			cur_xs <- t(cross$genotypes$real[rownames(cur_map)[which(cur_map[,1]==j)],])
			corM <- cor(cur_ys,cur_xs,use="pairwise.complete.obs")
			#cat(dim(cur_xs),"\n",dim(cur_ys),"\n",dim(result),"\n",dim(corM),"\n",dim(result[colnames(cur_ys),colnames(cur_xs)]),"\n")
			result[colnames(cur_ys),colnames(cur_xs)] <- corM
		}
	}
	invisible(result)
}

