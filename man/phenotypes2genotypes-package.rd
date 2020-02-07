\name{pheno2geno-package}
\alias{pheno2geno}
\docType{package}
\title{
	Pheno2Geno - High-throughput generation of genetic markers and maps from molecular phenotypes
}
\description{
Pheno2geno is an R package to generate genetic markers and maps out from molecular phenotypes. 
Currently supported breeding schemes are: Recombinant Inbred Lines (RIL), F2 and 
backcross (BC). 

The most important functions:
	\itemize{
    \item \code{\link{read.population}} - Reads the files into R.
		\item \code{\link{find.diff.expressed}} - Using Rank Product or student t-test analysis to select differentially expressed genes.
		\item \code{\link{scan.qtls}} - Scanning population data for qtls for use in cross.saturate function.
		\item \code{\link{generate.biomarkers}} - Converts continous gene expression measurments into discrete genetic markers.
    \item \code{\link{cross.denovo}} - Create de novo genetic map or vector showing how chromosomes should be assigned.
    \item \code{\link{cross.saturate}} - Saturate an existing genetic map with phenotype-derived markers.
  }
}
\details{
Background
Genetic markers and maps are instrumental for quantitative trait locus (QTL) mapping in segregating
populations. The resolution of QTL localisation depends on the number of informative recombinations in the
population and how well these recombinations are tagged by markers. Thus larger populations and denser
marker maps do a better job. Ideally there are enough markers to pinpoint all informative recombinations in the
population. In practice marker maps are often still too sparse. However, maps can be saturated or even be
derived de-novo from high-throughput gene expression, protein or metabolite abundance data. A fraction of
these molecular traits may show a clear multimodal distribution due to a major QTL effect and can therefore be
converted into useful genetic markers.
Results
We developed the pheno2geno R package for high-throughput generation of genetic markers and maps from
molecular phenotypes. Pheno2geno selects suitable phenotypes that show clear differential expression in the
1
founders. Mixture modelling is used to select phenotypes showing segregation ratios close to the expected
mendelian segregation ratios and transform these phenotypes into genetic markers suitable for map construction
and/or saturation. We demonstrate our method on 164 individuals from an A. thaliana Recombinant Inbred
Line (RIL) population. We show that pheno2geno is able to saturate the existing genetic map decreasing the
average distance between markers from 7.1 cM to 0.70 cM, pinpointing all recombinations in the population.
Using pheno2geno we created a de-novo map from the gene expression data that is twice as dense as the
original genetic map consisting of AFLP markers.
}

\author{ 
	Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{k.zych@rug.nl}
}
\references{
  
  Arends D., Zych K., Li Y., van der Velde K. J., Joosen R. V. L., Ligterink W. and Jansen R. C. (2013):
  Pheno2Geno - High-throughput generation of genetic markers and maps from molecular phenotypes. \bold{in press}
  
  Breitling, R.; Armengaud, P.; Amtmann, A.; and Herzyk, P.(2004) Rank Products:A simple, 
  yet powerful, new method to detect differentially regulated genes in replicated microarray experiments, \bold{FEBS Letter}, 57383-92

  Broman, K. W. and Sen,
  \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. (2009) \emph{A
  guide to QTL mapping with R/qtl.}  \bold{Springer}.  \url{http://www.rqtl.org/book}
}

\keyword{ package }

\seealso{
  \itemize{
    \item{\code{\link{read.population}}}{ - Load genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{find.diff.expressed}}}{ - Using Rank Product or student t-test analysis to select differentially expressed genes.}
    \item{\code{\link{scan.qtls}}}{ - Scanning population data for qtls for use in cross.saturate function.}
    \item{\code{\link{cross.denovo}}}{ - Create de novo genetic map or vector showing how chromosomes should be assigned.}
    \item{\code{\link{cross.saturate}}}{ - Saturate existing map.}
  }
}