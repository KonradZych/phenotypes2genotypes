\name{readFiles}
\alias{readFiles}

\title{Loading genotype and phenotype data}

\description{
  Loads genotype, phenotype, genetic map data files into R environment into a population object.
}

\usage{
	readFiles(offspring="offspring",founders="founders",map="maps",founders_groups,verbose=FALSE,debugMode=0)
}

\arguments{
 \item{offspring}{ Core used to specify names of children phenotypic ("core_phenotypes.txt") and genotypic ("core_genotypes.txt") files.}
 \item{founders}{ Core used to specify names of parental phenotypic ("core_phenotypes.txt") file. }
 \item{map}{ Core used to specify names of genetic ("map_genetic.txt") and physical ("map_physical.txt") map files. }
 \item{founders_groups}{ Specify groups of individuals in founders data, see description below and \code{\link[RankProd]{RP}} for more details }
 \item{verbose}{ Be verbose}
 \item{debugMode}{ 1: Print out checks, 2: print additional time information }
}

\value{
  Object of class population. See \code{\link{createPopulation}} for more details about structure.
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
Genotype file should be basically the same as phenotype file, but, apart from row/colnames, only three values are possible: 0,1 or NA. Example of genotype file structure:
  \tabular{lrrrrr}{
                      \tab "individual1"   \tab "individual2"   \tab "individual3"     \tab "individual4"     \tab "individual5"   \cr
"marker"                    \tab 1    \tab 1      \tab 0      \tab 1    \tab 0    \cr
"marker2"                   \tab NA      \tab 1    \tab 0      \tab 1      \tab 0    \cr
"marker3"             \tab 1      \tab 1    \tab 1      \tab 1      \tab 0    \cr
"marker4"        \tab 1      \tab NA    \tab 1      \tab 1      \tab 0    \cr
"marker5"               \tab NA    \tab 1    \tab 1    \tab 1    \tab 0    \cr
}
Map files should have really simple structure, always three columns, no header. First column contains rownames, second - chromosome number and third - position on chromosome (in cM for genetic or Mbp for physical map).
Secodn and third column can contain only numbers (any NA, Inf, etc, will cause dropping of file). Rownames should match either ones from genotype file or ones from phenotype file, depending which one you want to use 
map with (see toGenotypes for more information). Example of map file structure:
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
	Konrad Zych \email{konrad.zych@uj.edu.pl}, Danny Arends \email{Danny.Arends@gmail.com}
	Maintainer: Konrad Zych \email{konrad.zych@uj.edu.pl}
}

\examples{
	\dontrun{
	### simplest call possible
	population <- readFiles(founders_groups=c(0,0,0,1,1,1))
	### more informative one
	population <- readFiles(founders_groups=c(0,0,0,1,1,1),verbose=TRUE,debugMode=1)
	### imagine you prefer parents and children instead of founders and offspring:
	population <- readFiles(offspring="children",founders="parents",founders_groups=c(0,0,0,1,1,1)verbose=TRUE,debugMode=1)
	### etc.. when you load it, you may want to inspect it:
	population$founders$phenotypes[1:10,]
	}
}

\seealso{
  \code{\link{findDiffExpressed}}
  \code{\link{toGenotypes}}
  \code{\link{createPopulation}}
  \code{\link{intoPopulation}}
}

\keyword{manip}
