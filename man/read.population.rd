\name{read.population}
\alias{read.population}

\title{Loading genotype and phenotype data}

\description{
  Loads genotype, phenotype, genetic map data files into R environment into a population object.
}

\usage{
  read.population(offspring = "offspring", founders = "founders", map = "map", 
  foundersGroups, populationType = c("riself", "f2", "bc", "risib"), 
  readMode = c("normal","HT"), threshold=0.05, verbose = FALSE, debugMode = 0, 
  n.cluster=1, ...)
}

\arguments{
 \item{offspring}{ Core used to specify names of children phenotypic ("core_phenotypes.txt") genotypic ("core_genotypes.txt") and annotations ("core_annotations.txt") files.}
 \item{founders}{ Core used to specify names of parental phenotypic ("core_phenotypes.txt") file. }
 \item{map}{ Core used to specify names of genetic ("map_genetic.txt") and physical ("map_physical.txt") map files. }
 \item{foundersGroups}{ Specify groups of individuals in founders data, see description for more details }
   \item{populationType}{ Type of the population data was obtained from:
   \itemize{
    \item{riself}{ - RILs by selfing.}
    \item{f2}{ - f2 cross.}
    \item{bc}{ - back cross.}
    \item{risib}{ - RILs by sibling mating.}
  } 
 }
 \item{readMode}{HT, or High-Throughput mode should be used when the very large dataset is processed (at least 10000 probes). Then files are read in chunks intead of at once.
  To avoid R memory limits, only probes showing differential expression between parent are selected. Size of the chunk and threshold for assesing significance can be specified
  (see description of ... parameter).}
 \item{threshold}{ - threshold for assesing probes that are differentially expressed between parents. 0.05 by default.}
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
 \item{n.cluster}{ number of cores used for calcuations }
 \item{...}{ Parameters passed to high-throughtput function:    \itemize{
    \item{transformations}{ - how should the data be transformed (see \code{\link{transformation}})}
    \item{sliceSize}{ - number of lines to be read at once byt HT function. 5000 by default.}
  }}
}

\value{
 An object of class \code{\link{population}}.
}

\details{
  Function is working on tab delimited files. 
  Phenotype files, both for founders and offspring, should have header, containing column names (so names of individuals). All the other rows should start with rowname (unique).
  Rownames and colnames are only values allowed to be not numeric. After file is read into R, check is performed and rows and columns containing values that are not numeric and not convertable to numeric, will be removed
  from dataset. Rownames should match between founders and offspring. After loading founders file in, all non-matching rows are removed. Example of phenotype file structure:
  \tabular{lrrrrr}{
                      \tab "individual1"   \tab "individual2"   \tab "individual3"     \tab "individual4"     \tab "individual5"   \cr
"marker"                    \tab 8.84494695336781    \tab 9.06939381429179      \tab 9.06939381429179      \tab 7.72431126650435   \tab 6.04480152688572    \cr
"marker2"                   \tab 9.06939381429179      \tab 7.85859536346299    \tab 8.84494695336781      \tab 6.04480152688572      \tab 7.72431126650435    \cr
"marker3"             \tab 6.04480152688572      \tab 6.04480152688572    \tab 7.85859536346299      \tab 7.72431126650435      \tab 7.85859536346299    \cr
"marker4"        \tab 6.04480152688572     \tab 7.85859536346299    \tab 6.04480152688572      \tab 8.84494695336781      \tab 7.85859536346299    \cr
"marker5"               \tab 7.72431126650435    \tab 7.72431126650435    \tab 17.85859536346299    \tab 7.85859536346299   \tab 7.85859536346299    \cr
}
Genotype file should have basically the same structure as the phenotype file. The genotypes codes are exactly the same as in r/qtl - for F2 populations:
AA - 1, AB - 2, BB - 3, not BB - 4, not AA - 5, missing - NA and for BC and RILs: AA - 1, BB - 2, missing - NA (see \code{\link[qtl]{read.cross}} for details.) 
Example of genotype file structure:
  \tabular{lrrrrr}{
                      \tab "individual1"   \tab "individual2"   \tab "individual3"     \tab "individual4"     \tab "individual5"   \cr
"marker"                    \tab 1    \tab 1      \tab 2      \tab 1    \tab 2    \cr
"marker2"                   \tab NA      \tab 1    \tab 2      \tab 1      \tab 2    \cr
"marker3"             \tab 1      \tab 1    \tab 1      \tab 1      \tab 2    \cr
"marker4"        \tab 1      \tab NA    \tab 1      \tab 1      \tab 2    \cr
"marker5"               \tab NA    \tab 1    \tab 1    \tab 1    \tab 2    \cr
}
Map files should have really simple structure, always three columns, no header. First column contains rownames, second - chromosome number and third - position on chromosome (in cM for genetic or Mbp for physical map).
Secodn and third column can contain only numbers (any NA, Inf, etc, will cause dropping of file). Rownames should match either ones from genotype file or ones from phenotype file, depending which one you want to use 
map with (see generate.biomarkers for more information). Example of map file structure:
  \tabular{lrr}{
"marker"                    \tab 1    \tab 0       \cr
"marker2"                   \tab 1      \tab 1.2     \cr
"marker3"             \tab 1      \tab 1.2      \cr
"marker4"        \tab 1      \tab 2     \cr
"marker5"               \tab 1    \tab 3      \cr
}
You have also to specify groups ion founders file, so which columns come from which parent. Let's imagine, you have measured both parents in triplo and data for first parent is in columns 1,3 and 5, for second parent - columns 2,4,6.
Founders groups should be c(0,1,0,1,0,1) then. Always use only 0 and 1 to specify groups.
}

\author{
  Konrad Zych \email{k.zych@rug.nl}, Danny Arends \email{Danny.Arends@gmail.com}
  Maintainer: Konrad Zych \email{k.zych@rug.nl}
}

\examples{
  \dontrun{
  ### simplest call possible
  population <- read.population(founders_groups=c(0,0,0,1,1,1))
  ### more informative one
  population <- read.population(founders_groups=c(0,0,0,1,1,1),verbose=TRUE,debugMode=1)
  ### imagine you prefer parents and children instead of founders and offspring:
  population <- read.population(offspring="children",founders="parents",
    founders_groups=c(0,0,0,1,1,1),verbose=TRUE,debugMode=1)
  ### etc.. when you load it, you may want to inspect it:
  population$founders$phenotypes[1:10,]
  }
}

\seealso{
  \itemize{
    \item{\code{\link{add.to.population}}}{ - Adding data to existing population object.}
    \item{\code{\link{create.population}}}{ - Create new object of class population.}
  }
}

\keyword{manip}
