# demgsa: An R package for carrying out global sensitivity analysis of demographic models

`demgsa` is an R package that provides functions to carry out a global 
sensitivity analysis (GSA) of demographic models created via the 
[RAMAS GIS](http://ramas.com/software.htm) Software. `demgsa` includes 
functions to make easier each major step of a GSA:

1. Creation of new demographic models parameterized with random sets of values,
which represent the uncertainty and/or stochasticity inherent in each parameter.
2. Creation of scripts to batch run RAMAS GIS models.
3. Reading and collating of model inputs and results for further analysis.

Currently the `demgsa` package is available via M. Aiello-Lammens' 
[GitHub account](https://github.com/mlammens). The package repository is
[here](https://github.com/mlammens/demgsa),
and we encourage Clones, Forks, and Pull Requests.
The easiest way to install this package is to use the 
[devtools](https://github.com/hadley/devtools) package. However, the 
package vignette will not be created upon install. An html version
of the vignette, which functions as a `demgsa` tutorial, is 
[here](http://htmlpreview.github.io/?https://github.com/mlammens/demgsa/blob/master/inst/doc/demgsa_tutorial.html).

## Install demgsa via devtools::install_github

```
## Check for devtools. Install if not already installed.
if( !( "devtools" %in% installed.packages() ) ){
  install.packages( "devtools" ) 
}

## Load devtools package
library( devtools )

## Check for demgsa. Install if not already installed
if( !( "demgsa" %in% installed.packages() ) ){
  devtools::install_github( "mlammens/demgsa" ) 
}

## Load demgsa package
library( demgsa )
```


