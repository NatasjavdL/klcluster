# (k,l)-clustering
This package contains an R library for clustering of trajectories.
It has been developed as part of a Msc Thesis. The initial idea of (k,l)-clustering comes from the following paper: ["Approximating (k,ℓ)-center clustering for curves"](https://arxiv.org/abs/1805.01547)

The data attached to the package is public data from R. Mann. For more information see ["Landscape complexity influences route-memory formation in navigating pigeons"](http://rsbl.royalsocietypublishing.org/content/10/1/20130885).

For an introduction to the package I would like to refer you to the vignettes and R documentation.

This package is developed under the GLP 2.0 license in accordance with other packages that are required.


## Installing the package

To install the package, first please install the devtools package: `install.packages(devtools)`.
Next you can easily install this package by using the install_github function: `install_github("NatasjavdL/klcluster")`.

Note that it was built on a Windows machine. In case you are installing on Windows you might get an error:
```
Error: loading failed
Execution halted
*** arch - x64
ERROR: loading failed for 'i386'
```

Then try the following:
```
library(devtools)   
options(devtools.install.args = "--no-multiarch")   
install_github("NatasjavdL/klcluster")   
```
