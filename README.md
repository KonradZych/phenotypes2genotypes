The pheno2geno package
=================
R package to create genetic maps out of gene expression data. Currently supported breeding schemes are: Recombinant inbred line, F2 and backcross.
The master branch is the current stabile version and is tagged. The master branch can be installed into R, any development branches can contain errors.

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

However, then you will install current version, that could be under development. If you want the error and warning free version, click "Files to download", selected latest
version (currently 0.4.4), download and unpack zip file and then install it:

	 $ R CMD INSTALL phenotypes2genotypes
	


Starting
--------
Load the library in the R interface by the following command (in R):
    
    $ > library(pheno2geno)                                   # Load the library
    $ > ?pheno2geno                                           # Show the help

To read in data files:
    
    $ > setwd(paste(.Library,"pheno2geno/data",sep="/"))
    $ > population <- readFiles(founders_groups=c(0,0,1,1))             
    $ WARNING: There is no genotypic file for rils: children_genotypes.txt genotypic data for rils will be simulated
    $ WARNING: There is no map file for rils: children_probes.txt further processing will take place without taking it into account


To start RankProd analysis:

    $ > population <- findDiffExpressed(population)
	

Data files
--------
Our workflow is based on assumption that data files provided to our software have certain names (can be modified by user) and format (unmodifiable). 
To be sure your analysis runs smoothly please acknowledge data format specification included in readFiles help file in DETAILS section, you can access by:
    
	$ > ?readFiles


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

e) Check it:

    $ R CMD check phenotypes2genotypes

f) If it's warning-free, go back to b) or submit a patch.

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
