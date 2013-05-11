\name{create.population}
\alias{create.population}
\alias{population}

\title{Create a population object}

\description{
  Create a new population object from phenotype data already loaded in the R environment
}

\usage{
  create.population(offspringPhenotypes, founders, foundersGroups, offspringGenotypes, mapsGenetic, mapsPhysical,
  populationType=c("riself", "f2", "bc", "risib"), noWarn=FALSE, verbose=FALSE, debugMode=0)
}

\arguments{
  \item{offspringPhenotypes}{ A matrix that contains the phenotype data measured on the offspring (required).}
  \item{founders}{ A matrix that contains the phenotype data measured on the founders (required).}
  \item{foundersGroups}{ When multiple measurement for the founders are present this is used to group the founders. The format is a matrix that contains the phenotype data measured on the (required).}
  \item{offspringGenotypes}{ Matrix containing any known offspring genotype data (optional).}
  \item{mapsGenetic}{ Matrix containing a known genetic map (optional).}
  \item{mapsPhysical}{ Matrix containing a known physical map (optional).}
  \item{populationType}{ Type of population the expression data was obtained from:
    \itemize{
      \item{riself}{ - Recombinant inbred line by selfing.}
      \item{f2}{ - F2 cross.}
      \item{bc}{ - Back cross.}
      \item{risib}{ - Recombinant inbred line by sibling mating.}
    }
  }
  \item{noWarn}{ If TRUE, no warnings will be produced. }
  \item{verbose}{ Be verbose. }
  \item{debugMode}{ Either use 1 or 2, this will modify the amount of information returned to the user. 1) Print out checks, 2) Print additional time information.}
}

\value{
  An object of class \code{\link{population}}. 
  This is a complex object containing all the information needed for the entire pheno2geno analysis. It's structure looks like 
  below (depending on which optional information was supplied):
  \itemize{
    \item{$offspring}{ - Section in the object which contains all data related to the offspring:
      \itemize{
        \item{$phenotypes}{ - Offspring phenotype data, stored as a numeric matrix, Rows: phenotypes, Columns: individuals.}
        \item{$genotypes}{ - (Optional) Offspring genotype data:
          \itemize{
            \item{$real}{ - The original data when a known genetic map is provided by the user - numeric matrix, Rows: markers, Columns: individuals.}
            \item{$simulated}{ - Simulated genetic map produced by \code{\link{generate.biomarkers}} from phenotype data - numeric matrix, Rows: markers, Columns: individuals.}
          }
        }
      }
    }
    \item{$founders}{ - Section in the object which contains all data related to the founders:
      \itemize{
        \item{$phenotypes}{ - Founders phenotype data, stored as a numeric matrix, Rows: phenotypes, Columns: individuals.}
        \item{$groups}{ - Groups a founder belong to when replicates are available, storad as a vector of 0s and 1s, specifying per column which founder phenotype data belongs to which group.}
        \item{$RP}{ - Results from the t.test or RankProd analysis on the founders phenotype data, by \code{\link{find.diff.expressed}.}
      }
    }
    }
    \item{$maps}{ - Section in the object which contains all data related to maps:
      \itemize{
        \item{$genetic}{ Genetic map, stored as a numeric matrix, Rows: markers, 2 Columns [1] Chromosome, [2] Position on chromosome in cM.}
        \item{$physical}{ Physical map, stored as anumeric matrix, Rows: markers, 2 Columns [1] Chromosome , [2] Position on chromosome in Mbp.}
      }
    }
  }
}

\details{
  When all required information is provided an object of class \code{\link{population}} is returned.
}

\author{
	Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{k.zych@rug.nl}
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
    \item{\code{\link{add.to.population}}}{ - Adding data to an existing population object.}
  }
}

\keyword{manip}
