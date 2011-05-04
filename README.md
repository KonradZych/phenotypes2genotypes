The pheno2geno package
=================
This package provides you with possibility to create genetic map out of expression data for children and parents. Working currently for RILs, preparing to run Back Cross analysis.
Mind the fact that this is still pre-release version, so not everything is perfect. Also use rather versions that were tagged, you can install them and perform R check without warnings.

Dependencies
------------
R software environment from [www.r-project.org](http://www.r-project.org/ "www.r-project.org"), qtl package from [www.rqtl.org](http://www.rqtl.org, "www.rqtl.org") and RankProd 
package from [www.bioconductor.org](http://www.bioconductor.org/packages/release/bioc/html/RankProd.html "www.bioconductor.org")

Installation
------------
Prepare your environment by following these two steps:

- Download and Install the R environment
- Install the qtl and RankProd packages into the R environment

Then install into R by using (from a terminal / commandline):

    $ git clone git://github.com/KonradZych/phenotypes2genotypes.git  # Download the repository
    $ R CMD INSTALL phenotypes2genotypes                            # Install the package

However, then you wull install current version, that could be under development. If you want the error and warning free version, click "Files to download", selected latest
version (currently 0.4.4), download and unpack zip file and then install it:

	 $ R CMD INSTALL phenotypes2genotypes
	
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
	

Data files
--------
Our workflow is based on assumption that data files provided to our software have certain names (can be modified by user) and format (unmodifiable). 
To be sure your analysis runs smoothly please acknowledge data format specification included in readFiles help file in DETAILS section, you can access by:
    
	$ > ?readFile


TODO
--------------------
TODO


Contributing
------------

Want to contribute? Great!

a) Clone it:

    $ git clone git://github.com/KonradZych/phenotypes2genotypes.git 
	
b) Install it:

    $ R CMD INSTALL phenotypes2genotypes
	
c) Run it:

    $ > library(pheno2geno)
	
d) Modify some code. (Search -> 'TODO')

e) Go back to 2, or

f) Submit a patch

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
