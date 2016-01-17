## ------------------------------------------------------------------------
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

## ---- echo=FALSE---------------------------------------------------------
data( snpl_nocc_results )
data( snpl_2m_results )

## ------------------------------------------------------------------------
system.file("extdata", package = "demgsa")

## ---- eval=FALSE---------------------------------------------------------
#  ## Set working directory -- change this line according to where you placed
#  ## these files!
#  setwd( "inst/extdata/" )
#  
#  ## Call the sensitivity function using the sens_config file for nocc
#  sensitivity( "sens_config_nocc.txt" )

## ----eval=FALSE----------------------------------------------------------
#  ## Get a list of mp files to extract results from
#  mp_list <- list.files( pattern = "snpl_nocc.*mp", full.names = TRUE )
#  
#  ## Call mp.mult.results
#  mp.mult.results( mp.file.list = mp_list, out.csv = "snpl_nocc_results.csv", spatial = TRUE )

## ---- eval=FALSE---------------------------------------------------------
#  ## Read in the resulting CSV file
#  snpl_nocc_results <- read.csv( "snpl_nocc_results.csv" )

## ------------------------------------------------------------------------
summary( snpl_nocc_results )

## ------------------------------------------------------------------------
plot( data = snpl_nocc_results, exp.min.n ~ GrowthRt )

## ---- eval=FALSE---------------------------------------------------------
#  ## Call the sensitivity function using the sens_config file for 2m slr
#  sensitivity( "sens_config_2m.txt" )

## ----eval=FALSE----------------------------------------------------------
#  ## Get a list of mp files to extract results from
#  mp_list_2m <- list.files( pattern = "snpl_2m.*mp", full.names = TRUE )
#  
#  ## Call mp.mult.results
#  mp.mult.results( mp.file.list = mp_list_2m, out.csv = "snpl_2m_results.csv", spatial = TRUE )

## ---- eval=FALSE---------------------------------------------------------
#  ## Read in the resulting CSV file
#  snpl_2m_results <- read.csv( "snpl_2m_results.csv" )

## ------------------------------------------------------------------------
all( snpl_nocc_results$GrowthRt == snpl_2m_results$GrowthRt )
all( snpl_nocc_results$stdev.avg == snpl_2m_results$stdev.avg )

## ------------------------------------------------------------------------
plot( x = snpl_nocc_results$exp.min.n, y = snpl_2m_results$exp.min.n,
      xlab = "EMA for NoCC scenario", ylab = "EMA for 2m SLR scenario" )
abline( a = 0, b = 1, col="red" )

