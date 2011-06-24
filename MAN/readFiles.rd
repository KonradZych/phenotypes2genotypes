\name{readFiles}
\alias{readFiles}

\title{Reading geno/phenotypic files.}

\description{
  Reads geno/phenotypic/map files into R environment into special object.
}

\usage{
	readFiles(offspring="offspring",founders="founders",map="maps",sep="",verbose=FALSE,debugMode=0)
}

\arguments{
 \item{offspring}{ Core used to specify names of children phenotypic ("core_phenotypes.txt") and genotypic ("core_genotypes.txt") files.}
 \item{founders}{ Core used to specify names of parental phenotypic ("core_phenotypes.txt") file. }
 \item{map}{ Core used to specify names of genetic ("map_genetic.txt") and physical ("map_physical.txt") map files. }
 \item{sep}{ Separator of values in files. Passed directly to read.table, so "" is a wildcard meaning whitespace.}
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
}

\value{
  Object of class population. See createPopulation for more details about structure.
}

\details{
  Function is working on tab delimited files. First row contains colnames, first column - rownames. No remarks, headers etc. If you fish to use files of
  your own format, please check createPopulation and intoPopulation functions.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
	Under tender patronage of: Danny Arends \email{Danny.Arends@gmail.com}
}

\examples{
	\dontrun{
	setwd(paste(.Library,"pheno2geno/data",sep="/"))
	population <- readFiles()
	population$founders$phenotypes[1:10,]
	}
}

\seealso{
  \code{\link{preprocessData}}
  \code{\link{toGenotypes}}
  \code{\link{createPopulation}}
  \code{\link{intoPopulation}}
}

\keyword{manip}
