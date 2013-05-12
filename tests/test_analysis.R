require(pheno2geno)
children <- read.csv(file="offspring_phenotypes.csv",header=TRUE,row.names=1)
parents <- read.csv(file="parental_phenotypes.csv",header=TRUE,row.names=1)
genotypes <- read.csv(file="genotypes.csv",header=TRUE,row.names=1)
map <- read.csv(file="map.csv",header=TRUE,row.names=1)

#with parental data
population <- create.population(children,parents,c(0,0,0,0,0,0,1,1,1,1,1,1),genotypes,verbose=TRUE)
population <- find.diff.expressed(population)
population <- generate.biomarkers(population, threshold=0.1, verbose=T, debug=2)
population <- scan.qtls(population,verbose=T,step=0, map="physical")

####THREE WAYS TO ASSIGN CHROMOSOMES
set.seed(101010)
cross_newmap <- cross.denovo(population,n.chr=16,map="physical",comparisonMethod=sumMajorityCorrelation,reOrder=TRUE,use.orderMarkers=FALSE,verbose=TRUE,debugMode=2)
cross_newmap <- smooth.geno(cross_newmap,3)

cross_saturated <- cross.saturate(population,map="physical",verbose=TRUE,debugMode=2)
cross_saturated <- smooth.geno(cross_saturated,3)

#To output folder
write.cross(cross_newmap,file="yeast_cross_newmap.csv",format="csv")
write.cross(cross_saturated,file="yeast_cross_enriched.csv",format="csv")


#Enriched versus Old
png("enriched_vs_old.png",width=1200,height=1200)
markersCorPlot(cross_saturated,population,"physical",comparisonMethod=sumMajorityCorrelation)
dev.off()

#New versus Old
png("new_vs_old.png",width=1200,height=1200)
markersCorPlot(cross_newmap,population,"physical",comparisonMethod=sumMajorityCorrelation)
dev.off()

#Enriched versus Old comparison
png("enriched_vs_old_comparison.png",width=1200,height=1200)
plotMapComparison(cross_saturated,population,"physical")
dev.off()
