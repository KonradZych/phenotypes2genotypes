\name{pheno2geno-package}
\alias{pheno2geno}
\docType{package}
\title{
	Pheno2geno - generating genetic markers and maps from molecular phenotypes.
}
\description{
Pheno2geno is an R package to generate genetic markers and maps out of phenotype expression data. 
Currently supported breeding schemes are: Recombinant Inbred Lines (RIL), F2 and 
backcross (BC). 

The most important functions:
	\itemize{
        \item \code{\link{read.population}} - Reads geno/phenotypic files into R.
		\item \code{\link{find.diff.expressed}} - Uses Student t-test to select differentially expressed genes.
		\item \code{\link{generate.biomarkers}} - Converts continous gene expression measurments into discrete genetic markers.
		}
}
\details{
TODO
}
\author{ 
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}
\references{
  
  Zych, K.; Arend, D. and Jansen R.C. (2011) Pheno2geno: R package for creation of genetic maps out of gene expression data, \bold{in press}
  
  Breitling, R.; Armengaud, P.; Amtmann, A.; and Herzyk, P.(2004) Rank Products:A simple, 
  yet powerful, new method to detect differentially regulated genes in replicated microarray experiments, \bold{FEBS Letter}, 57383-92

  Broman, K. W. and Sen,
  \if{latex}{\out{\'S}}\if{html}{\out{&#346;}}\if{text}{S}. (2009) \emph{A
  guide to QTL mapping with R/qtl.}  \bold{Springer}.  \url{http://www.rqtl.org/book}
}

\keyword{ package, qtl }

\seealso{
TODO
}