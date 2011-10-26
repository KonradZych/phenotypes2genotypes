setwd("D:/data/yeast_new/step1")
require(pheno2geno)
bremcross <- read.cross("csvr",file="yeast_brem_cross.csv",geno=c(0,1))
bremcross <- convert2riself(bremcross)
children <- read.csv(file="offspring_phenotypes.csv",header=TRUE,row.names=1)
parents <- read.csv(file="parental_phenotypes.csv",header=TRUE,row.names=1)
genotypes <- read.csv(file="genotypes.csv",header=TRUE,row.names=1)
map <- read.csv(file="map.csv",header=TRUE,row.names=1)

population <- createPopulation(children,parents,c(0,0,0,0,0,0,1,1,1,1,1,1),t(genotypes),maps_physical=map)
population <- findDiffExpressed(population)
population <- toGenotypes(population, treshold=0.01,verbose=T,debug=2)
set.seed(101010)

####THREE WAYS TO ASSIGN CHROMOSOMES
#cross_newmap <- assignChromosomes(cross_newmap,population,"physical",how=correlationRule,reOrder=TRUE,verbose=TRUE)
#cross_newmap <- assignChromosomes(cross_newmap,population,"physical",how=majorityRule,reOrder=TRUE,verbose=TRUE)
cross_newmap <- assignChromosomes(population,80,"physical",how=sumMajorityRule,reOrder=TRUE,verbose=TRUE)
cross_newmap <- smoothGeno(cross_newmap,3)

cross_enriched <- enrichExistingMap(population,map="physical",verbose=TRUE,debugMode=2)
cross_enriched <- smoothGeno(cross_enriched,3)


#To output folder
write.cross(cross_newmap,file="yeast_cross_newmap.csv",format="csv")
write.cross(cross_enriched,file="yeast_cross_enriched.csv",format="csv")


#Enriched versus Old
png("enriched_vs_old.png",width=1200,height=1200)
markersCorPlot(cross_enriched,population,"physical",show=mean)
dev.off()

#New versus Old
png("new_vs_old.png",width=1200,height=1200)
markersCorPlot(cross_newmap,population,"physical",show=mean)
dev.off()

#Enriched versus Old comparison
png("enriched_vs_old_comparison.png",width=1200,height=1200)
plotMapComparison(cross_enrich,population_a,"physical")
dev.off()
