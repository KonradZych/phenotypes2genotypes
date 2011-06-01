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
	intoRil(ril=NULL, dataObject, dataType=c("parental","children"), selectedRows=NULL, selectedCols=NULL)
}

\arguments{
 \item{ril}{ object of class ril, data should be put into, by default - new object will be created}
 \item{dataObject}{ matrix of data to be put into ril object}
 \item{dataType}{ what kind of data dataObject contains (parental phenotypic, children phenotypic)}
 \item{selectedRows}{ rows to be selected from dataObject, by default - all}
 \item{selectedCols}{ cols to be selected from dataObject, by default - all }
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
	ril <- intoRil(dataObject=parental, dataType="parental")
	ril <- intoRil(ril, children, "children")
}

\seealso{
  \code{\link{readFiles}}
}

\keyword{manip}
