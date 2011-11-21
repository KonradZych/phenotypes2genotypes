\name{findDiffExpressed}
\alias{findDiffExpressed}
\alias{groupLabels}

\title{Finding differentially expressed genes.}

\description{
  Using Rank Product analysis to select differentially expressed genes.
}

\usage{
	findDiffExpressed(population,use=c("ttest","rankprod"),verbose=FALSE,debugMode=0,...)
}

\arguments{
 \item{population}{ an object of class \code{\link{population}}}
 \item{use}{ which method should be used for selecting differentially expressed probes:
  \itemize{
    \item{ttest}{ - student t-test}
    \item{rankprod}{ - Rank Product using \code{\link{RP}} function}
  } 
  }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
 \item{...}{ Additional arguments passed to RP function. }
}

\value{
  Object of class \code{\link{population}}.
}

\details{
  TODO
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\references{
  Hong F, Breitling R, McEntee CW, Wittner BS, Nemhauser JL, Chory J. RankProd: a bioconductor package for detecting differentially expressed genes in meta-analysis. Bioinformatics. 2006 Nov 15;22(22):2825-7.
}

\examples{
	data(yeastPopulation)
	yeastPopulation <- findDiffExpressed(yeastPopulation)

}

\seealso{
	\itemize{
    \item \code{\link[RankProd]{RP}} - Add description
    \item \code{\link{readFiles}} - Add description
    \item \code{\link{toGenotypes}} - Creating genotypes from children phenotypes
    \item \code{\link{showRPpval}} - Printing out p-values calculated by the findDiffExpressed function
    \item \code{\link{plotRPpval}} - Plotting p-values calculated by the findDiffExpressed function
  }
}

\keyword{manip}
