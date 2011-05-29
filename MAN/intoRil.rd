\name{intoRil}
\alias{intoRil}
\alias{parental}
\alias{parentalRows}
\alias{parentalCols}
\alias{children}
\alias{childrenRows}
\alias{childrenCols}

\title{intoRil}

\description{
  Putting data into ril object.
}

\usage{
	intoRil(ril=NULL, parental=NULL, parentalRows=NULL, parentalCols=NULL, children=NULL, childrenRows=NULL, childrenCols=NULL)
}

\arguments{
 \item{ril}{ object of class ril, data should be put into, by default - new object will be created}
 \item{parental}{ matrix of parental data to be put into ril object}
 \item{parentalRows}{ rows to be selected from parental, by default - all}
 \item{parentalCols}{ cols to be selected from parental, by default - all }
 \item{children}{ matrix of parental data to be put into ril object}
 \item{childrenRows}{ rows to be selected from children, by default - all }
 \item{childrenCols}{ cols to be selected from children, by default - all }
}

\value{
  Object of class ril.
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
	parental <- read.table("parental_phenotypes.txt",sep="")
	children <- read.table("children_phenotypes.txt",sep="")
	ril <- intoRil(parental=parental, children=children)
}

\seealso{
  \code{\link{readFiles}}
}

\keyword{manip}
