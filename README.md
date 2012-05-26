The pheno2geno package
======================
Pheno2geno is an R package to create genetic maps out of phenotype expression data. 
Currently supported breeding schemes are: Recombinant Inbred Lines (RIL), F2 and 
backcross (BC). The master branch is the current stable version and is tagged. The 
master branch can be installed into R, development branches may contain errors.

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

```
    $ git clone git://github.com/KonradZych/phenotypes2genotypes.git  # Download the repository
    $ R CMD INSTALL phenotypes2genotypes                            # Install the package
```

However, then you will install current version, that could be under development. If you want the error and warning free version, click "Files to download", selected latest
version (currently 0.4.4), download and unpack zip file and then install it:

```
    $ R CMD INSTALL phenotypes2genotypes
```

Starting
--------
Load the library in the R interface by the following command (in R):

```R
    $ > library(pheno2geno)                                   # Load the library
    $ > ?pheno2geno                                           # Show the help
```

To read in data files, use the readFiles function:

```R    
    $ > population <- readFiles(founders_groups=c(0,0,1,1))
```

We start by finding differentialy expressed phenotypes, by using wither RankProd 
or a T-test:

```R
    $ > population <- findDiffExpressed(population)
```

Data files
--------
Our workflow is based on assumption that data files provided to our software have 
certain names (can be modified by user) and format (unmodifiable). To be sure your 
analysis runs smoothly please acknowledge data format specification included in 
readFiles help file in DETAILS section, you can access by:

```R
    $ > ?readFiles
```

TODO
----

See inst/TODO.txt

Contributing
------------

Want to contribute? Great!

a) Clone it:

```
    $ git clone git://github.com/KonradZych/phenotypes2genotypes.git 
```

b) Install it:

```
    $ R CMD INSTALL phenotypes2genotypes
```

c) Run it:

```R
    $ > library(pheno2geno)
```

d) Modify some code. (Search -> 'TODO')

e) Check it:

```
    $ R CMD check phenotypes2genotypes
```

f) If it's warning-free, go back to b) or submit a patch.

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
