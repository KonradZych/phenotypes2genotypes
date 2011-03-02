\name{phenotypes2genotypes-package}
\alias{phenotypes2genotypes}
\alias{ph2g}
\docType{package}
\title{
	Basic QTL - very simple qtl mapping and genotype markers ordering of messed data
}
\description{
	Includes the following functionality:
	\itemize{
    \item \code{\link{qtlAnalysis}} - All in one, just give it wd, file names and output name for plot image
    \itemize{
	  \item \code{\link{heatmapqtl}} - Basic QTL mapping currently done using qtlbyttest for every marker in dataset
      \item \code{\link{qtlbyttest}} - Single marker mapping by using t-test statistics
      \item \code{\link{pathway}} - Basic pathway creation using single marker
    }
    \item \code{\link{un_neighbor}} - Ordering markers from messed genotypical data
	\itemize{
      \item \code{\link{un_drop_markers}} - Removing non-informational markers, that are higly corelated with more than specified percentage of others
    }
    \item \code{\link{utilities}} - Ploting routines and helper functions
	\itemize{
      \item \code{\link{persp_qtl_map}} - produces nice and highly informative perspective plot of QTL map
	  \item \code{\link{makebinary}} - Produces binary matrix from dataset by splitting it by specified treshold
	  \item \code{\link{un_recombination}} - Calculates recombination factor between two markers
    }
	}
}
\details{
TODO
}
\author{ 
	Konrad Zych \email{konrad.zych@uj.edu.pl}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
	Under tender patronage of: Danny Arends \email{Danny.Arends@gmail.com}
}
\references{
TODO
}
\keyword{ package,qtl }

\seealso{
TODO
}