\name{pheno2geno-package}
\alias{pheno2geno}
\docType{package}
\title{
	Tools for the construction of genetic maps from phentype data
}
\description{
	Includes the following functionality:
	\itemize{
        \item \code{\link{read.population}} - Reading geno/phenotypoc files into R.
		\item \code{\link{find.diff.expressed}} - Using Rank Prod to select differentially expressed genes.
		\item \code{\link{generate.biomarkers}} - Convert a phenotypematrix into suitable genotypes and save it into R/qtl cross object.
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