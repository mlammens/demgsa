sensitivity <- function( sens.config.file ) {
# RAMAS Metapopulation Sensitivity Analysis Program
# Author: Matthew Aiello-Lammens
# Date: 28 January 2011
# Update: 29 August 2011; 3 November 2011, 3 December 2011, 8 December 2011
#
# Args:
#  sens.config.file: File that contains the configureations for a sensitivity analysis run, based on 
#  the inlcuded sens.config.txt file.
#
# Returns:
#
#  No direct returns to the R consule
#
# This script is written to read two files that are associated with a single
# species.  One file is assumed to be the 'high' estimates of parameter values
# and the other file is assumed to be the 'low' estimates of the parameter
# values.
#
# NOTES: In this version we are interested in examing the effects of varying the following
# parameters:
# - survival rates
# - fecundity rates
# - survival variability
# - fecundity variability
# - dispersal rates
# - inter-population correlation
# - Rmax (for contest density dependence)
# Density dependence will be selected by the user, and it is recommended that
# the user examine variability for all feasible dens. dep. scenarios (i.e., 
# perform a sensitivity analysis for both ceiling and contest type dens. dep.)
###################################################################################################
# Defined functions used in script
# ----------------------------------------------------------------------------------------------- #
# matrix.check( mat1, mat2) - checks that element-for-element, the members of mat1 are less than
#  those for mat2
matrix.check <- function( mat1, mat2 ){
  if ( all( mat2 >= mat1 ) ) { return(TRUE) }
  else { return(FALSE) }
}
# ----------------------------------------------------------------------------------------------- #
# rand.st.matr( corrMat, low.matr, high.matr, const.matr)
# Create a new stage (or standard deviation) matrix based on a low and high matrix and the user defined 
# inter-dependence structure (correlation) among and within fecundity and survival parameters
#
# Args:
#  corrMatt: A matrix defining the inter-dependence structutre among and within fecundity and survival parameters; 
#    E.g., 1, 1, 0, where 1 is completely dependent and 0 is completely independent. This parameter is defined
#    in the sens.config.txt file.  More information is availble in the default sens.config.txt file.
#
#  low.matr: Stage (st.dev.) matrix of all low estimates
#
#  high.matr: Stage (st.dev.) matrix of all high estimates
#
#  const.matrix: Constraints matrix - associated with RAMAS Metapop.  Indicates which stage matrix elements
#    are fecundity parameters and which are survival parameters
#
# Returns:
#  StMatNew: A new stage (st.dev.) matrix that has been varied between the low and high values based on the inter-dependence structure
#
rand.st.matr <- function( corrMat, low.matr, high.matr, const.matr, st.rvs ) {
  survCor <- corrMat[1]
  fecCor <- corrMat[2]
  sxfCor <- corrMat[3]

  # Number of stages should be equivalent to the number of rows or columns of the stage matrix. This
  # is not correct if there is sex structure.
  stages <- nrow(low.matr)
  
  # Determine the number of RandUnifs and repeats for survival 
  if ( survCor == 0 ) {
    surv.runifs <- stages
    surv.reps <- 1
  } else if ( survCor == 1 ) {
    surv.runifs <- 1
    surv.reps <- stages
  } else {
    stop( 'Interdependence value of variablity for survival is not equal to 0 or 1.\n This program is not yet equiped for other values.')
  }
  # Determine the number of RandUnifs and repeats for fecundity
  if ( fecCor == 0 ) {
    fec.runifs <- stages
    fec.reps <- 1
  } else if ( fecCor == 1 ) {
    fec.runifs <- 1
    fec.reps <- stages
  } else {
    stop( 'Interdependence value of variablity for fecundity is not equal to 0 or 1. This program is not yet equiped for other values.')
  }

  if ( sxfCor == 0 ){
    # Vary survival and fecundity seperately
    
    # Fecundity
    randFec <- rep( st.rvs[1:fec.runifs],fec.reps )
    st.rvs <- st.rvs[-(1:fec.runifs)]
    randFecMat <- matrix( rep( randFec, stages ), nrow=stages, byrow = TRUE )
    StMatNewFec <- (( high.matr - low.matr ) * randFecMat ) + low.matr 
    # Eliminate non-fecundity stage elements, as noted by the constraints matrix
    StMatNewFec <- StMatNewFec * (1 - const.matr)
    
    # Survival
    randSurv <- rep( st.rvs[1:surv.runifs],surv.reps ); 
    st.rvs <- st.rvs[-(1:fec.runifs)]
    randSurvMat <- matrix( rep( randSurv, stages ), nrow=stages, byrow = TRUE )
    StMatNewSurv <- (( high.matr - low.matr ) * randSurvMat ) + low.matr
    # Eliminate non-survival stage elements, as noted by the constraints matrix
    StMatNewSurv <- StMatNewSurv * const.matr
    
    # Add Fec and Surv matrices to make new matrix
    StMatNew <- StMatNewFec + StMatNewSurv
    
  } else if ( sxfCor == 1){
    # Vary survival and fecundity dependently
    randMat <- matrix( rep( st.rvs[1:surv.runifs],(stages * surv.reps) ), nrow=stages, byrow=TRUE )
    StMatNew <- (( high.matr - low.matr ) * randMat ) + low.matr    
  }
  return( StMatNew )
}
###################################################################################################
# Start main program
# ******************
#
# ----------------------------------------------------------------------------------------------- #
# Read values from the sensitivity analysis config file (e.g., sens.config.txt)
# ----------------------------------------------------------------------------------------------- #
# Read config file if it exists
if( file.exists( sens.config.file ) ) {
  sens.config <- readLines( sens.config.file )
} else {
  stop( paste("Could not find sensitivity configuration file: ",sens.config.file) )
}
# Check that file is not empty, and if it passes check, make new variable without comment lines
if( length( sens.config ) > 0 ) {
  # Make new variable without comment lines
  sens.config.vals <- sens.config[ -(grep("#",sens.config)) ]
  # Remove blank spaces
  sens.config.vals <- gsub( " ", "", sens.config.vals)
} else {
  stop( paste("Sensitivity configuration file is empty: ", sens.config.file ) )
}
# Check that the correct number of input parameters is in the configuration file
# PROGRAMING NOTE: req.in.vals is a variable that may need to be changed or adjusted if additional
# inputs are written into the sensitivity.r script.
req.in.vals <- 13
if ( length( sens.config.vals ) < req.in.vals ) {
  stop( paste("Missing input parameters in sensitivity configuration file: ", sens.config.file ) )
}
# 
# Read mp.file.names - A vector of *.mp files used in the sensitivity analysis.  The first element of the vector
# MUST be an *.mp file that contains all of the low estimates of parameter values and the last element 
# of the vector MUST be an *.mp file that conatins all of the high estimates of parameter values.
mp.file.names.line <- grep('mp.file.names', sens.config.vals)
mp.file.names <- unlist(strsplit( sens.config.vals[ mp.file.names.line ], split='=' ))
mp.file.names <- unlist(strsplit( mp.file.names[2], split=',' ))
#
# Read surv.fec.corr - A matrix defining the inter-dependence structutre among and within fecundity and 
# survival parameters; E.g., 1, 1, 0 where 1 is completely dependent and 0 is completely 
# independent
surv.fec.corr.line <- grep('surv.fec.corr', sens.config.vals)
surv.fec.corr <- unlist(strsplit( sens.config.vals[ surv.fec.corr.line ], split='=' ))
surv.fec.corr <- as.numeric(unlist(strsplit(surv.fec.corr[2], split=',' )))
# Check that surv.fec.corr has 3 elements
if ( length( surv.fec.corr ) != 3 ) {
  stop( "Incorrect survival-fecundity interdependence structure (surv.fec.corr)" )
}
#
# Read sens.iter - The desired number of iterations of the sensitibity analysis - this program
# will produce sens.iter number of new *.mp files
sens.iter.line <- grep('sens.iter', sens.config.vals)
sens.iter <- unlist(strsplit( sens.config.vals[ sens.iter.line ], split='=' ))
sens.iter <- as.numeric( sens.iter[2] )
#
# Read out.dir - The base directory for where the new *.mp files should go
out.dir.line <- grep('out.dir', sens.config.vals)
out.dir <- unlist(strsplit( sens.config.vals[ out.dir.line ], split='=' ))
out.dir <- out.dir[2]
# 
# Read out.name - The base name for the new *.mp files
out.name.line <- grep('out.name', sens.config.vals)
out.name <- unlist(strsplit( sens.config.vals[ out.name.line ], split='=' ))
out.name <- out.name[2]
#
# Read batch.file: Base name for newly created batch file(s)
batch.file.line <- grep('batch.file', sens.config.vals)
batch.file <- unlist(strsplit( sens.config.vals[ batch.file.line ], split='=' ))
batch.file <- batch.file[2]
#
# Read bat.file.cnt: The number of batch files to be created.  Using this feature allows the user
# to create multiple batch files, which when executed, will run on multiple processors if available
bat.file.cnt.line <- grep('bat.file.cnt', sens.config.vals)
bat.file.cnt <- unlist(strsplit( sens.config.vals[ bat.file.cnt.line ], split='=' ))
bat.file.cnt <- as.numeric( bat.file.cnt[2] )
#
# Read rand.samp: 
rand.samp.line <- grep('rand.samp', sens.config.vals)
rand.samp <- unlist(strsplit( sens.config.vals[ rand.samp.line ], split='=' ))
rand.samp <- rand.samp[2]
#
# Read rv.file: Name of the file in which either random variables will be saved to or from which random
# variables will be taken from
rv.file.line <- grep('rv.file', sens.config.vals)
rv.file <- unlist(strsplit( sens.config.vals[ rv.file.line ], split='=' ))
rv.file <- rv.file[2]
#
# Read use.rv.file: TRUE/FALSE Use random variables from rv.file, rather than writing to the rv.file
use.rv.file.line <- grep('use.rv.file', sens.config.vals)
use.rv.file <- unlist(strsplit( sens.config.vals[ use.rv.file.line ], split='=' ))
use.rv.file <- as.logical( use.rv.file[2] )
#
# Read use.sad: TRUE/FALSE Use stable age (stage) distribution when calculated the stage initial abundance values
use.sad.line <- grep('use.sad', sens.config.vals)
use.sad <- unlist(strsplit( sens.config.vals[ use.sad.line ], split='=' ))
use.sad <- as.logical( use.sad[2] )
# 
# Read pop.disp.depend: 0 or 1, A value indicating whether uncertainty of population dispersal values
# is completely independent (0) or dependent (1) among populations
pop.disp.depend.line <- grep('pop.disp.depend', sens.config.vals)
pop.disp.depend <- unlist(strsplit( sens.config.vals[ pop.disp.depend.line ], split='=' ))
pop.disp.depend <- as.numeric( pop.disp.depend[2] )
#
# Read pop.disp.dch.include: TRUE/FALSE - Include DCH files when varying dispersal scenarios.  See sens.config.txt for more
# information regarding this option
pop.disp.dch.include.line <- grep('pop.disp.dch.include', sens.config.vals)
pop.disp.dch.include <- unlist(strsplit( sens.config.vals[ pop.disp.dch.include.line ], split='=' ))
pop.disp.dch.include <- as.logical( pop.disp.dch.include[2] )
#
# Read pop.kch.include: TRUE/FALSE - Include kch files when varying dispersal scenarios.  See sens.config.txt for more
# information regarding this option
pop.kch.include.line <- grep('pop.kch.include', sens.config.vals)
pop.kch.include <- unlist(strsplit( sens.config.vals[ pop.kch.include.line ], split='=' ))
pop.kch.include <- as.logical( pop.kch.include[2] )
#
# Read user.defined.dd: TRUE/FALSE - Include user defined d-d parameters stored in PopData_df.
# information regarding this option
user.defined.dd.line <- grep('user.defined.dd', sens.config.vals)
user.defined.dd <- unlist(strsplit( sens.config.vals[ user.defined.dd.line ], split='=' ))
user.defined.dd <- as.logical( user.defined.dd[2] )
#
# ----------------------------------------------------------------------------------------------- #
# Use mp.read.r function to read *.mp files and convert into a sorted 'list' structure.
# For more information on the sorted mp.read list structure, see the 'mp.read.r' function
# file.
# First check to see that at least two *.mp files are given by the user.
if( length( mp.file.names) < 2 ){
  stop("User Input Error: User must provide at least two *.mp file names in input value mp.file.names")
}
# mp.files will be a nested list structure.  At the top level, each list element will be a seperate
# *.mp file converted to a list structure by the mp.read function.
mp.files <- lapply( mp.file.names, mp.read )

# For simplicity of coding and to work with existing code, set the first and last elements of the
# mp.files structure as mp.low and mp.high, respectively
mp.low <- mp.files[[1]]
mp.high <- mp.files[[ length(mp.files) ]]
# Determine Metapop version number for each file and check for match
mp.lowVersion <- mp.low$version
mp.highVersion <- mp.high$version
if ( mp.lowVersion != mp.highVersion ) {
  stop("Warning: The Metapop Versions do not match between these two files.")
}
# Assigning sorted list structure of the *.mp files to the variable names mp.low and mp.high
mp.low <- mp.low$mp.file
mp.high <- mp.high$mp.file

# ----------------------------------------------------------------------------------------------- #
# Set Variables that may be used multiple times in program
pop.num <- length(mp.low$PopList) # Determine the number of populations in the low mp file 
pop.num.high <- length(mp.high$PopList) # Number of populations in high mp file

# Create batch file path and base name
batch.file <- paste( out.dir, '/', batch.file, sep="")

# Check to see of batch files with batch.file in their name exist already
if( length( dir( path = out.dir, pattern = batch.file ))) {
  print( paste(batch.file, ' files may already exists! New lines well be appended to existing files.', sep="") )
}

# ----------------------------------------------------------------------------------------------- #
# Identify difference between mp files
# ************************************
# General flow to identify which elements of the mp lists do not match:
# - Create a vector of logical vals called 'noChange' - assuming low and high files are completely different 
#   (all noChange values set to FALSE). 
# - For each element of the mp lists, check if there is "no change" between the two lists 
# - Note the absence of changes as a logical value, setting that element in the noChange vector
#  as TRUE, otherwise leave the value as FALSE
# ----------------------------------------------------------------------------------------------- #

noChange <- logical(51) # Note: All values are initiated as FALSE
# Loop through all of the none 'list' elements of the mp 'list'. The 'list' elements are not 
# checked here because the == operator does not allow for logical testing of nested lists

# The number of 'list' type elements to check is different depending on whether
# the population numbers of the low and high mp files agree.  If they do not agree,
# then parameters relating to populations are ignored (i.e. not checked for differences)
if ( pop.num == pop.num.high ) {
  # non.list elements in the case that population numbers match
  non.list <- c(1:27,29:35,37,39:45,47:51)
} else {
  # non.list elements in the case that population numbers DO NOT match
  non.list <- c(1:27,29:30,32:33,35,37,40:44,47:51)
}
for (i in non.list ) {
  if ( all(mp.low[[i]] == mp.high[[i]]) ) {
    noChange[i] <- TRUE
  } 
}
# ----------------------------------------------------------------------------------------------- #
# Define a list of parameters that are not permited to vary (i.e., Non-variable Parameters), 
# this list is based on discussions with R. Ackakaya.
# Version 1.x Non-variable parameters: Number of stages in stage-matrix, Number of replications, 
#  Duration, Time-step size, Time-step units, Initial number of time-steps to exclude from risk calc
#  
# If some non-variable parameters do not agree then stop the sensitivity analysis
nonVariable <- c( 1, 2, 4, 19, 20, 51)
if ( all( noChange[nonVariable] ) ) {
  print("Non-Variable parameters agree!")
} else {
  print("ERROR: Some Non-Variable parameters do not agree between *.mp files. Non-variable parameters are: ")
  print( guide[nonVariable] )
  stop()
}
# ----------------------------------------------------------------------------------------------- #
# If number of populations is equal, check the PopList element: Check if ANY differences exist, if not, 
# then we can mark 'noChange[28]' as true and move on.  If false, we'll return to this later.
if ( pop.num == pop.num.high ) {
  if ( identical( mp.low$PopData_df, mp.high$PopData_df ) ) {
    noChange[28] <- TRUE
  }
} else {
  print('WARNING: The number of populations in mp.low and mp.high do not match. Ignoring population differences and using all mp.low population specific parameters.')
}
# ----------------------------------------------------------------------------------------------- #
# Check Stage-Matrix Differences. StMatr is element 36 in mp list structure:
#  - If there are multiple stage matrices, only the first one is used in this analysis.  Like wise
#  for st.dev. matrices
#   - Ignores differences in StMatrSurvMult and StMatrFecMult, as these values are used only to
#  create new matrices, and not actually used during a simulation run
#
# Check if there is more than one stage-matrix
StMatLowCnt <- length(mp.low$StMatr)
StMatHighCnt <- length(mp.high$StMatr)
if ( StMatLowCnt > 1 ) {
  print("Warning: mp.file.low contains greater than one stage matrix.  Only the first matrix will be used")
}
if ( StMatHighCnt > 1 ) {
  print("Warning: mp.file.high contains greater than one stage matrix.  Only the first matrix will be used")
}
# Considering only the first matrix, check if the two matrices are identical
StMatLow <- mp.low$StMatr[[1]]$Matr
StMatHigh <- mp.high$StMatr[[1]]$Matr
if ( all( StMatLow == StMatHigh ) ) {
  noChange[36] <- TRUE
}
# ----------------------------------------------------------------------------------------------- #
# Check Standard Deviation-Matrix Differences. 
# Check if there is more than one stdev matrix
SDMatLowCnt <- length(mp.low$SDMatr)
SDMatHighCnt <- length(mp.high$SDMatr)
if ( SDMatLowCnt > 1 ) {
  print("Warning: mp.file.low contains greater than one st.dev. matrix.  Only the first matrix will be used")
}
if ( SDMatHighCnt > 1 ) {
  print("Warning: mp.file.high contains greater than one st.dev. matrix.  Only the first matrix will be used")
}
# Considering only the first matrix, check if the two matrices are identical
SDMatrLow <- mp.low$SDMatr[[1]]$Matr; 
SDMatrHigh <- mp.high$SDMatr[[1]]$Matr;
if ( all( SDMatrLow == SDMatrHigh ) ) {
  noChange[38] <- TRUE
} 
#-------------------------------------------------------------------------------------------------- #
# Check that Stage Specific Properties (e.g., Name, Weight, Exclude, Basis for DD, Breeding Prop.)
# are that same.  In this current version, no variance is accounted for in these values.
StPropDiffs <- unlist( mp.low$StProp ) == unlist( mp.high$StProp )
if ( all(StPropDiffs) ) {
  noChange[46] <- TRUE
} else {
  print( "Warning: Stage specific properties, such as name, weight, exclude, basis for DD, etc. do not match!" )
}
#-------------------------------------------------------------------------------------------------- #

if ( all( noChange == TRUE ) ) {
  print("No differences between *.mp files were found (that this program looks for).  No sensitivity analysis will be performed.")
} else {
  # Begin creation of new *mp files for use in Sensitivity Analysis
  # ***************************************************************
  print( paste("Begin new *.mp file(s) creation. ", sens.iter, " files will be created.") ) 
  # Tell the user which parameters differ between the high and low *.mp files
  parDiffs <- which(noChange == FALSE)
  print( "Parameters that vary between *.mp files: " )
  print( guide[parDiffs] )

  # Determine the total number of random variables that will be needed in this analysis
  #   Initiate a list of 51 zeros
  potential.var.cnt <- rep(0,51)
  #   Now set the ones that are fixed
  #    Stage Initial Abundance 
  potential.var.cnt[45] <- 1
  
  # Determine Number Vars for Stage Matrix
  surv.cor <- surv.fec.corr[1]
  fec.cor <- surv.fec.corr[2]
  sxf.cor <- surv.fec.corr[3]
  # Determine the number of vars for survival
  if ( surv.cor == 0 ) {
    surv.vars <- mp.low$Stages
  } else if ( surv.cor == 1 ) {
    surv.vars <- 1
  }
  # Determine the number of vars for fecundity
  if ( fec.cor == 0 ) {
    fec.vars <- mp.low$Stages
  } else if ( fec.cor == 1 ) {
    fec.vars <- 1
  }
  # Determine the total number of variables for the stage matrix variation
  if ( sxf.cor == 1 ) {
    if ( surv.cor == 1 ) {
      stage.vars <- 1
    } else {
      stage.vars <- surv.vars
    }
  } else if ( sxf.cor == 0 ) {
    stage.vars <- surv.vars + fec.vars
  }
  print( paste('Count of *independently* varying stage parameters:', stage.vars) )
  potential.var.cnt[36] <- stage.vars # The number of rvs for stage matrix
  
  # Assume that since the stage matrix and st.dev. matrix currently share the same correlation structure
  # that we can use the stage.vars as the st.dev.vars
  st.dev.vars <- stage.vars
  potential.var.cnt[38] <- st.dev.vars # Number of rvs for st.dev. matrix

  # Calculate the number of additional variables if population numbers match
  if ( pop.num == pop.num.high ) {
    potential.var.cnt[34] <- 1 # If the correlation matrix is used, only one rv is used
    potential.var.cnt[33] <- 1 # If correlation-distance function is used

    # Dispersal
    # 
    # If using DCH files, then vars = 1, otherwise check pop.disp.depend
    if ( pop.disp.dch.include ) {
      potential.var.cnt[31] <- 1
    } else { 
      if ( pop.disp.depend == 1) {
        potential.var.cnt[31] <- 1 # Dispersal varied using on rv value
      } else if ( pop.disp.depend == 0 ) {
        potential.var.cnt[31] <- pop.num*pop.num # Dispersal varied using unique rv value for each element of disp matrix
      } else {
        stop( 'Inappropriate value for pop.disp.depend used.  Value must be 0 or 1.')
      } # End 'if pop.disp.depend' conditional
    } # End 'if pop.disp.dch.includ' conditional
    
    # Stage initial abundance values
    # If the user selects use.sad, then there are no rvs added, otherwise one is added
    if( !use.sad ) { 
      potential.var.cnt[45] <- 1 
    }
    
    # Population specific information:
    # - One for initial abundance 
    # - One for MaxR
    # - One for User defined density dependent parameters
    potential.var.cnt[28] <- 3
    # Add one additional if pop.kch.include is true
    if( pop.kch.include ){ potential.var.cnt[28] <- 3 }
  } # End 'if pop.num == pop.num.high' loop
  
  # Tally up the total number of parameters that are potetnailly being varied
  var.cnt <- sum( potential.var.cnt )
  print( paste( 'Count of *independently* varying *.mp file parameters (including stage parameters):', var.cnt ))
  ###print('DEBUG OUTPUT: sens.iter and var.cnt'); print(sens.iter); print(var.cnt) ### DEBUG LINE
  
  if ( use.rv.file ) {
    # Import random variables from the rv.file
    uni.rv.mat <- read.csv( rv.file )
    uni.rv.mat <- as.matrix( uni.rv.mat )
    print('Debug output: uni.rv.mat'); print(uni.rv.mat) ### DEBUG LINE
    print( dim(uni.rv.mat) )
    
    # Check if imported uni.rv.mat is the correct size
    if ( ! (all ( dim(uni.rv.mat) == c(sens.iter, var.cnt) ) ) ) {
      stop("Random variables from rv.file do not match the dimensions required for this sensitivity run.")
    }  
  } else if ( !use.rv.file ) {
    # Create new random variables and save to a new rv.file
    
    # Sample all of the uniform random variables based on value of 'rand.samp'
    # that will be needed using the Latin Hypercube 
    # Sampling
    if ( rand.samp == 'lhs' ) {
      # Use Latin Hypercube Sampling
      uni.rv.mat <- randomLHS( sens.iter, var.cnt )
    } else if ( rand.samp == 'urand' ) {
      # Use uniform random variables without any constraints
      uni.rv.mat <- matrix( runif( sens.iter*var.cnt ), sens.iter, var.cnt )
    } else if ( rand.samp == 'custom' ) {
      # Create a custom random variable sampling
      uni.rv.mat <- matrix( rep( .5, (sens.iter*var.cnt) ), sens.iter, var.cnt )
    } else {
      stop("Incorrect random sampling scheme chosen for input parameter 'rand.var'.")
    }
    
    # Save random variables to rv.file
    write.csv( uni.rv.mat, file=rv.file, row.names=FALSE )    
  } # End 'if use.rv.file'
  
  #----------------------------------------------------------------------------------------------#
  # Begin 'for' loop for creation of *.mp files. 'sens.iter' number of loops will be carried out
  for ( iter in 1:sens.iter ) {
    print('***************************************************************')
    print( paste( 'Sensitivity iteration:', iter ) ) 
    # Define a new mp.read list structure with the same values as the mp.low structure
    mp.new <- mp.low
    # Get this iterations uniform random variables. As these variables are used, remove them from uni.rv.
    uni.rv <- uni.rv.mat[iter,]
    print('Random numbers used to vary *.mp parameters:'); print(uni.rv) 
    #----------------------------------------------------------------------------------------------#  
    # Stage Matrix (Fecundity and Survival) 
    # *************************************
    if ( noChange[36] == FALSE ) {
      # Generate a new stage matrix
      if ( !matrix.check( StMatLow, StMatHigh ) ){
        print("WARNING: Elements of the stage matrix in the 'Low' file are greater than those in the 'High' file.")
      }
      # Call 'rand.st.matr' function to create a new stage matrix
      mp.new$StMatr[[1]]$Matr <- rand.st.matr( surv.fec.corr, StMatLow, StMatHigh, mp.new$ConstraintsMatr, uni.rv[1:stage.vars] )
      uni.rv <- uni.rv[-(1:stage.vars)]
      print('New Stage Matrix'); print(mp.new$StMatr[[1]]$Matr)       
    } else {
    print("No differences between stage matrices were observed.  No new stage matrix created.")
    }
    ###print('uni.rv after stg.matrix creation'); print(uni.rv) ### DEBUG LINE
    #----------------------------------------------------------------------------------------------#  
    # Standard Deviation Matrix (Fecundity and Survival) 
    # **************************************************
    if ( noChange[38] == FALSE ) {
      # Generate a new standard deviation matrix
      if ( !matrix.check( SDMatrLow, SDMatrHigh ) ){
        print("WARNING: Elements of the standard deviation matrix in the 'Low' file are greater than those in the 'High' file.")
      }
      # Call 'rand.st.matr' to create new std.dev matrix, using the same inter-dependence structure as fecundity and survival
      mp.new$SDMatr[[1]]$Matr <- rand.st.matr( surv.fec.corr, SDMatrLow, SDMatrHigh, mp.new$ConstraintsMatr, uni.rv[1:stage.vars] )
      print('New Std. Dev. Matrix'); print(mp.new$SDMatr[[1]]$Matr)
      uni.rv <- uni.rv[-(1:stage.vars)]
      ###print('uni.rv after st.dev matrix creation'); print(uni.rv) ### DEBUG LINE
    } else {
    print("No differences between standard deviation matrices were observed.  No new std. dev. matrix created.")
    }
    #----------------------------------------------------------------------------------------------#  
    # If the number of populations match between the low and high *.mp files, then vary the population
    # specific values, including the dispersal matrix, correlation matrix, pop specific information, 
    # and stage intial abundance values.
    if ( pop.num == pop.num.high ) {
      #----------------------------------------------------------------------------------------------#  
      # Population Specific Information
      # *******************************
      if ( noChange[28] == FALSE ) {
        # Copy two data.frame(s) of all population specific information to new variables
        PopLow_df <- mp.low$PopData_df
        PopHigh_df <- mp.high$PopData_df
        # Set new mp file data frame as low data frame
        PopNew_df <- PopLow_df
        
        # Change initial population abundance (N0).  All populations varied using same
        # random uniform number draw.
        InitAbund_Low <- PopLow_df$InitAbund
        InitAbund_High <- PopHigh_df$InitAbund
        InitAbund_New <- ((InitAbund_High - InitAbund_Low)*uni.rv[1]) + InitAbund_Low
        PopNew_df$InitAbund <- as.integer(InitAbund_New)
        uni.rv <- uni.rv[-1]
        ####print('uni.rv after pop specific InitAbund'); print(uni.rv)
        
        # Change Rmax. All populations varied using the same random uniform number draw.
        MaxR_Low <- PopLow_df$MaxR
        MaxR_High <- PopHigh_df$MaxR
        MaxR_New <- ((MaxR_High - MaxR_Low)*uni.rv[1]) + MaxR_Low
        PopNew_df$MaxR <- MaxR_New
        uni.rv <- uni.rv[-1]
        ####print('uni.rv after maxR'); print(uni.rv)
        
        ## -------------------------------------------------------------------- ##
        # User defined density dependence parameters
        ###browser()
        if ( user.defined.dd ){
          # Determine UDD pars
          udd_pars <- grep(pattern='udd',names(PopLow_df))
          udd_pars_new <- ((PopHigh_df[udd_pars]-PopLow_df[udd_pars])*uni.rv[1]) + PopLow_df[udd_pars]
          PopNew_df[udd_pars] <- udd_pars_new
          uni.rv <- uni.rv[-1]
        } else {
          ### FRAL specific addition ###
          # Remove the random variable generated to change UDD pars
          uni.rv <- uni.rv[-1]
        }
        
        # Set mp.new$PopList to PopNew_df
        # Note: this changes the structure of this list element, from a nested list of lists to a 
        #  data frame object
        mp.new$PopList <- PopNew_df
      } else {
      # Create a data.frame for the population data for the mp.new file using just the mp.low parameters if no differences were found.
        print( "No differences in population specific parameters were found.  Using mp.low values for the new *.mp file." )
        mp.new$PopList <- mp.low$PopData_df
      }
      #----------------------------------------------------------------------------------------------#
      # Dispersal Matrix
      # ****************
      if ( noChange[31] == FALSE | noChange[30] == FALSE ) {
        # First check to see if the mp files agree on UseDispDistFucn, if they do agree
        # and there is no inclusion of dch files, then only change the dispersal distance
        # function parameters. NOTE: if inclusion of the dch files in the sensitivity 
        # analysis is desired, then a dispersal distance function must be created.
        if ( mp.low$UseDispDistFunc & mp.high$UseDispDistFunc & !pop.disp.dch.include ) {
          mp.new$UseDispDistFunc <- TRUE
          DispDistFunc_New <- (( mp.high$DispDistFunc - mp.low$DispDistFunc) * uni.rv[1] ) + mp.low$DispDistFunc
          mp.new$DispDistFunc <- DispDistFunc_New
          print('New dispersal distance function parameters:'); print(mp.new$DispDistFunc)
          uni.rv <- uni.rv[-1]
        } else if ( pop.disp.dch.include ) { 
          # If we want to vary the dch file scenarios, then complete the following:
          # Select one of the dispersal scenarios from the list of input *.mp files

          # Set UseDispDistFunc as 'FALSE' in the new *.mp file.
          mp.new$UseDispDistFunc <- FALSE
          
          # Get the random value already drawn for changing dispersal. This value is between 0 and 1
          disp.matr.rand <- uni.rv[1]
          # Convert to an integer between 1 and the total number of mp files given by user
          disp.matr.num <- ceiling( disp.matr.rand * length(mp.files) )
          print( paste('Dispersal Matrix From MP File: ', mp.file.names[ disp.matr.num ] ) )
          
          # Set new dispersal matrix as the dispersal matrix of the randomly selected dispersal matrix
          mp.new$DispMatr <- mp.files[[ disp.matr.num ]]$mp.file$DispMatr
          # Set new dispersal distance function to the randomly selected mp file
          mp.new$DispDistFunc <- mp.files[[ disp.matr.num ]]$mp.file$DispDistFunc
          # Change associated DCH files to match those used in the dispersal scenario for the 
          # disp.matr.num selected.  Do this in both PopData_df and PopList
          mp.new$PopData_df$RelDisp <- mp.files[[ disp.matr.num ]]$mp.file$PopData_df$RelDisp
          mp.new$PopList$RelDisp <- mp.files[[ disp.matr.num ]]$mp.file$PopData_df$RelDisp
          uni.rv <- uni.rv[-1]
        } else if ( !all(mp.low$UseDispDistFunc,mp.high$UseDispDistFunc) & !pop.disp.dch.include ){
          # If one of the files includes a dispersal distance function AND we are NOT
          # interested in varying dch files, complete the following:

          # Set UseDispDistFunc as 'FALSE' in the new *.mp file.
          mp.new$UseDispDistFunc <- FALSE
          
          # Check matrix elements
          if ( !matrix.check( mp.low$DispMatr, mp.high$DispMatr ) ){
            print("WARNING: Elemnts of the dispersal matrix in the 'Low' file are greater than those in the 'High' file.")
          }
          # Generate new dispersal matrix, accounting for inter-dependence in uncertainty accross populations
          if ( pop.disp.depend == 1 ) {
            mp.new$DispMatr <- ((mp.high$DispMatr - mp.low$DispMatr)*uni.rv[1]) + mp.low$DispMatr
            uni.rv <- uni.rv[-1]
          } else {
            disp.rand.matr <- matrix( uni.rv[1:(pop.num*pop.num)], pop.num, pop.num )
            mp.new$DispMatr <- ((mp.high$DispMatr - mp.low$DispMatr)*disp.rand.matr) + mp.low$DispMatr
            uni.rv <- uni.rv[-(1:(pop.num*pop.num))]
          }
        } else {
          stop( "ERROR in parsing dispersal distance measures. Check dispersal related settings in config file." )
        } # End dispersal scenario 'if' conditional
      } else {
        print("No differences in the 'low' and 'high' dispersal matrices were observed. No new disp. matrix created.")
      } # End Dispersal Matrix 'if' loop
      #----------------------------------------------------------------------------------------------#  
      # Correlation Matrix
      # ******************
      if ( noChange[34] == FALSE | noChange[33] == FALSE) {
        # Check to see if either file uses a user-defined correlation matrix.  If not,
        # vary the parameters in the user defined matrix.  If so, vary the elements of the
        # matrix using a single uniform random variable.
        if ( mp.low$UseCorrDistFunc & mp.high$UseCorrDistFunc ) {
          mp.new$UseCorrDistFunc <- TRUE
          CorrDistFunc_New <- ((mp.high$CorrDistFunc - mp.low$CorrDistFunc) * uni.rv[1]) + mp.low$CorrDistFunc
          mp.new$CorrDistFunc <- CorrDistFunc_New
          print('New correlation distance function parameters:'); print( mp.new$CorrDistFunc )
          uni.rv <- uni.rv[-1]
          ####print('uni.rv after correlation function'); print(uni.rv)
        } else {
          mp.new$UseCorrDistFunc <- FALSE
          CorrLow <- mp.low$CorrMatr
          CorrHigh <- mp.high$CorrMatr
          if ( matrix.check( CorrLow, CorrHigh ) ){
            # Generate a new *symetric* correlation matrix.
            # Only using one random variable! This means we are assuming inter-dependence between populations
            mp.new$CorrMatr <- ((CorrHigh - CorrLow)*uni.rv[1]) + CorrLow
            uni.rv <- uni.rv[-1]
            ####print('uni.rv after correlation matrix'); print(uni.rv)
          } else {
            stop("Elements of the correlation matrix in the 'Low' file are greater than those in the 'High' file.")
          } # End matrix.check 'if' loop
        } # End UseCorrDistFunc 'if' loop
      } else {
        print("No differences in the 'low' and 'high' correlation matrices were observed. No new corr. matrix created.")
      } # End Correlation Matrix 'if' loop
      #----------------------------------------------------------------------------------------------#  
      # Stage Initial Abundance values
      # ******************************
      # If the user specified 'use stable age distribution', then change the values of the first column
      # of the StInit matrix to -1.
      # If the user specified that these values should vary, then allow these values to vary interdependently
      if ( use.sad ) {
        # Set first column to -1, which tells RAMAS to use the Stable Age Distn
        mp.new$StInit <- mp.low$StInit
        mp.new$StInit[,1] <- -1
      } else {
        # Vary the Stage Initial Abundance values
        mp.new$StInit <- round( ((mp.high$StInit - mp.low$StInit)*uni.rv[1]) + mp.low$StInit )
        uni.rv <- uni.rv[-1]
        ####print('uni.rv after stage initial abundance'); print(uni.rv)
      } # End use.sad 'if' loop
      #----------------------------------------------------------------------------------------------#
      # Population Carrying Capacity Through Time
      # *****************************************
      # If pop.kch.include = TRUE, then chose the KCH files from one of the scenarios
      if ( pop.kch.include ) {
        # Get random variable
        kch.rand <- uni.rv[1]
        # Convert to an integer between 1 and the total number of mp files given by user
        kch.num <- ceiling( kch.rand * length(mp.files) )
        print( paste('KCH Files From MP File: ', mp.file.names[ kch.num ] ) )
        # Set kch files in new *.mp file to those from selected scenario
        mp.new$PopData_df$KchangeSt <- mp.files[[ kch.num ]]$mp.file$PopData_df$KchangeSt
        # Set kch file in PopList information
        mp.new$PopList$KchangeSt <- mp.new$PopData_df$KchangeSt
        
        uni.rv <- uni.rv[-1]
      }
      
      
      #----------------------------------------------------------------------------------------------#
      
      
    } # End 'if (pop.num == pop.num.high)' conditional

    #----------------------------------------------------------------------------------------------#  
    # Write the new *.mp file
    # ***********************
    # Create new file name
    new.out.name <- paste( out.name, iter, '.mp', sep = "" )
    mp.new.file <- paste( out.dir, '/', new.out.name, sep="")
    
    # Check that out.dir exist, create if not
    if ( !file.exists( out.dir ) ) {
      print( paste('Creating new output directory: ', out.dir, sep="") )
      dir.create( out.dir )
    } else {
      print( paste('Writing new *.mp files to output directory: ', out.dir, sep="") )
    }

    # Check if new *.mp file exists - if yes, then do not overwrite
    if ( file.exists( mp.new.file ) ) {
      print( paste('WARNING: ',mp.new.file,' already exists! Will not overwrite', sep="") )
    } else {
      print( paste('Writing new *.mp file: ', mp.new.file, sep="") )
      mp.write( mp.new, mp.lowVersion, mp.new.file )
      # If the file is being created on a *nix system, then convert
      # the file to Windows format (CRLF)
      if(.Platform$OS.type == "unix") {
        print( 'Converting *.mp file to Windows format' )
        system( paste(sens.base.dir,'nix2win.sh ',mp.new.file, sep="") )
      }
    }
        
    # Create new batch file
    bat.file.num <- iter %% bat.file.cnt
    b.file <- paste( batch.file, bat.file.num, '.bat', sep = "" )
    
    # Write a line to the Windows Batch script that will run all of the new metapop files
    batch.line <- paste( 'CALL RunMP', new.out.name, '/RUN=YES', sep = " " )
    ### FRAL specific line ###
    #batch.line <- paste( 'CALL RunLabMP', new.out.name, '/RUN=YES', sep = " " )
    print( paste( "Writing a new line to batch file: ", b.file ) )
    write( batch.line, file = b.file, append = TRUE )
    
    # Create new shell script file (for running on Linux Boxes)
    sh.file <- paste( batch.file, bat.file.num, '.sh', sep = "" )
    
    # Write a line to the Shell script that will run the new metapop files
    shell.line <- paste( 'wine Metapop.exe', new.out.name, '/RUN=YES', sep = " " )
    write( shell.line, file = sh.file, append = TRUE )
    Sys.chmod( sh.file, mode="0777" )

    ####print('last look at uni.rv in sens.iter for loop'); print(uni.rv) ### DEBUG LINE
  } # End sens.iter for loop
} # End if noChange condition
} # End sensitivity function
