The pheno2geno package
======================
Pheno2geno is an R package to create genetic maps out of phenotype expression data. 
Currently supported breeding schemes are: Recombinant Inbred Lines (RIL), F2 and 
backcross (BC). The master branch is the current stable version and is tagged. The 
master branch can be installed into R, development branches may contain errors.

Dependencies
------------
+ R software environment from [www.r-project.org](http://www.r-project.org/ "www.r-project.org")
+ qtl package from [www.rqtl.org](http://www.rqtl.org "www.rqtl.org") also available on [CRAN](http://cran.r-project.org/web/packages/qtl/index.html "http://cran.r-project.org/web/packages/qtl/index.html")
+ mixtools package available on [CRAN](http://cran.r-project.org/web/packages/mixtools/index.html "http://cran.r-project.org/web/packages/mixtools/index.html")

Optional:
RankProd package from [www.bioconductor.org](http://www.bioconductor.org/packages/release/bioc/html/RankProd.html "www.bioconductor.org")

Installation
------------
The easiest way is to use package devtools. Type in R:

    install.packages("devtools")
    library(devtools)
    install_github(repo="phenotypes2genotypes",username="KonradZych")
    library("pheno2geno")

Starting
--------
Load the library in the R interface by the following command (in R):

    library(pheno2geno)

You can always access the help files of the package or for any function by typing:

    ?pheno2geno
    ?function.name

Or:

    help(pheno2geno)
    help(function.name)

To read in data files, use the read.population function:

    population <- read.population(founders_groups=c(0,0,1,1))

In the help file of this function there is a description of the expected file formats. For more information see the manual.

TODO
----

See inst/TODO.txt

Contributing
------------

Want to contribute? Great!

a) Clone a local version of the Github repository to your own hard disk:

    git clone git://github.com/KonradZych/phenotypes2genotypes.git

b) Install it from the commandline by using the following command:

    R CMD INSTALL phenotypes2genotypes

c) Then start R and load the library to make the functions available:

    library(pheno2geno)

d) Modify some code. (Search -> 'TODO')

e) To check if the package is able to install in R:

    R CMD check phenotypes2genotypes

f) If it's warning-free, go back to b) or submit a push request!

You can also just post comments on code / commits.

Disclaimer
----------
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License,
version 3, as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the GNU
General Public License, version 3, for more details.

A copy of the GNU General Public License, version 3, is available
at [http://www.r-project.org/Licenses/GPL-3](http://www.r-project.org/Licenses/GPL-3 "GPL-3 Licence")

Copyright (c) 2010-2012 GBIC - Danny Arends, Konrad Zych, Ritsert C. Jansen
