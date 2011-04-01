\name{plotChildrenExpression}
\alias{plotChildrenExpression}
\alias{markers}

\title{Plotting routine for children expression data.}

\description{
  Plots children expression data (boxplot) max value of parental exptression (red trangle) min (blue triangle) and mean(line) for selected markers.
}

\usage{
	plotChildrenExpression(ril, markers=1:100)
}

\arguments{
 \item{ril}{ Ril type object, must contain parental phenotypic data.}
 \item{markers}{ Numbers of markers to be plotted. }
}

\value{
	Plot.
}

\details{
  TODO
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
	Under tender patronage of: Danny Arends \email{Danny.Arends@gmail.com}
}

\examples{
	setwd(paste(.Library,"pheno2geno/data",sep="/"))
	ril <- readFiles()
	plotChildrenExpression(ril)
}

\seealso{
  \code{\link{plotChildrenExpression}}
}

\keyword{manip}
