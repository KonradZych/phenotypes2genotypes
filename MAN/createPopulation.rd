\name{create.population}
\alias{create.population}
\alias{population}

\title{Create a new population object}

\description{
  Create new population object from phenotype data
}

\usage{
	create.population(offspring_phenotypes, founders, founders_groups, offspring_genotypes, maps_genetic, maps_physical, populationType=c("riself", "f2", "bc", "risib"), no.warn=FALSE, verbose=FALSE,debugMode=0)
}

\arguments{
  \item{offspring_phenotypes}{ Matrix containing offspring phenotype data (have to be supported, if not - function quits with error).}
  \item{founders}{ Matrix containing founders phenotype data (have to be supported, if not - function quits with error).}
  \item{founders_groups}{ Specify groups im founders data (have to be supported, if not - function quits with error).}
 \item{offspring_genotypes}{ Matrix containing offspring genotype data (optional).}
 \item{maps_genetic}{ Matrix containing genetic map (optional).}
 \item{maps_physical}{ Matrix containing physical map (optional).}
  \item{populationType}{ Type of the population data was obtained from:
   \itemize{
    \item{riself}{ - RILs by selfing.}
    \item{f2}{ - f2 cross.}
    \item{bc}{ - back cross.}
    \item{risib}{ - RILs by sibling mating.}
  } 
 }
 \item{no.warn}{ If TRUE, no warnings will be produced.}
 \item{verbose}{ Be verbose.}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information.}
}

\value{
  An object of class \code{\link{population}}. 
  This is a complex object containing all the information needen for the pheno2geno analysis. It is splitted into few parts:
  \itemize{
    \item{$offspring}{ - containing all the offspring data:
      \itemize{
        \item{$phenotypes}{ - offspring gene expression (phenotype) data - numeric matrix, rows - markers, cols - individuals.}
        \item{$genotypes}{ - offspring genotype data:
          \itemize{
            \item{$real}{ - original data provided by the user - numeric matrix, rows - markers, cols - individuals.}
            \item{$simulated}{ - simulated by \code{\link{generate.biomarkers}} using phenotype data - numeric matrix, rows - markers, cols - individuals.}
          }
        }
      }
    }
    \item{$founders}{ - containing all the founders data:
      \itemize{
        \item{$phenotypes}{ - founders gene expression (phenotype) data - numeric matrix, rows - markers, cols - individuals.}
        \item{$groups}{ - vector of 0s and 1s, specifying which column in founders phenotype data belongs to which group.}
        \item{$RP}{ - results of t.test or RankProd analysis of founders phenotype data made by \code{\link{find.diff.expressed}.}
      }
    }
    }
    \item{$maps}{ - containing maps:
      \itemize{
        \item{$genetic}{ - genetic map - numeric array, rows - markers, 2 columns - 1 - chromosome marker lies on, 2 - position on chromosome in cM.}
        \item{$physical}{ - physical map - numeric array, rows - markers, 2 columns - 1 - chromosome marker lies on, 2 - position on chromosome in Mbp.}
      }
    }
  }
}

\details{
  An object of class \code{\link{population}}.
}

\author{
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	### simulating data
	population <- fake.population()
	offspring <- population$offspring$phenotypes
	founders <- population$founders$phenotypes
	founders_groups <- population$founders$groups
	population <- create.population(offspring,founders,founders_groups)
}

\seealso{
  \itemize{
    \item{\code{\link{read.population}}}{ - Loads genotype, phenotype, genetic map data files into R environment into a population object.}
    \item{\code{\link{add.to.population}}}{ - Adding data to existing population object.}
  }
}

\keyword{manip}
