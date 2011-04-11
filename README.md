The GoogleFinancial package
=================
This package provides a wrapper to download financial data from [Google Finance](http://www.google.com/finance "Google Finance"). 
It provides tools to model the downloaded financial data. This package provides automated loading of datasets 
from [Google Finance](http://www.google.com/finance "Google Finance"). It aims to predict stock prices 
for the upcoming 2 years . Furthermore the packages allows for exploration of the downloaded data like 
discovery of batch / localized effects. To provide easy data downloading and managment a pre-made list 
of companies from AEX, DowJones and NasDaQ is provided. However it is always possible to define your 
own set of companies / indices to monitor and study.

Dependencies
------------
R software environment from [www.r-project.org](http://www.r-project.org/ "www.r-project.org"), qtl package from [www.rqtl.org](http://www.rqtl.org, "www.rqtl.org") and RankProd 
package from [www.bioconductor.org](http://www.bioconductor.org/packages/release/bioc/html/RankProd.html "www.bioconductor.org")

Installation
------------
Prepare your environment by following these two steps:

- Download and Install the R environment
- Install the qtl and RankProd packages into the R environment

Then install into R by using (from a terminal / commanline):

    $ git clone git://github.com/KonradZych/phenotypes2genotypes.git  # Download the repository
    $ R CMD INSTALL phenotypes2genotypes                            # Install the package

Optionally you can install the pre-build packages by downloading the appropriate 
package for your operating system. 

Starting
--------
Load the library in the R interface by the following command (in R):
    
    $ > library(pheno2geno)                                   # Load the library
    $ > ?pheno2geno                                           # Show the help

To read in data files:
    
    $ > setwd(paste(.Library,"pheno2geno/data",sep="/"))
    $ > ril <- readFiles()             
    $ WARNING: There is no genotypic file for rils: children_genotypes.txt genotypic data for rils will be simulated
    $ WARNING: There is no map file for rils: children_probes.txt further processing will take place without taking it into account


To start RankProd analysis:

    $ > ril <- preprocessData(ril)

TODO
--------------------
TODO


Contributing
------------

Want to contribute? Great!

1. Clone it.
2. Compile it.
3. Run it.
4. Modify some code. (Search -> 'TODO')
5. Go back to 2, or
6. Submit a patch

You can also just post comments on code / commits.

Konrad Zych

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
Copyright (c) 2010 Danny Arends
