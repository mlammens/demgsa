mp.read <- function(mpFile) {
# Read data from the *.mp file.  Based on Metapop version 5 and 5.1 formats.
#
# Args:
#  mpFile: the name and path of the *.mp file to be read.  Path can be either from the 
#  current working directory or from the root directory
#
# Returns:
#  A list of two elements: version and mp.file
#
#  version: the Metapop version used to create the mp.file
#
#  mp.file: a sorted 'list' object containing all of the mp.file input information.  It does not 
#  contain any results information. Elements of the list strucutre are named to conform to naming 
#  scheme used in Metapop code.
#
###################################################################################################
# mp.file 'list' structure elements. Most element names from RAMAS Metapop source code
# The assignment of each element of the mp.file list is based on the Metapopulation input file
# type for version 5.0.  
#  1  MaxRep: maximum number of replications (integer; 1-10000)
#  2  MaxDur: Duration of simulation (integer; 1-500)
#  3  Demog_Stoch: Use demographic stochasticity? (boolean)
#  4  Stages: Number of stages (integer; 1-50)
#  5  Stages.IgnoreConstraints: Ingore constraints? (boolean - TRUE/FALSE)
#  6  Cat1: Information associated with catastrophe 1
#    This information is subdivided in RAMAS Input, but will not be subdivided
#    in the first version of the sensitivity analysis package
#  7  Cat2: Information associated with catastrophe 2
#    See note associated with Cat1
#  8  DDActing: Information regarding density dependence
#  9  Distrib: Distribution (Normal or Lognormal) to use for Environmental
#    Stochasticity.
#  10 AdvancedStochSetting: Advanced Environmental Stoch. Settings: (0 or PoolVar or 
#    NegSurvCorr)
#  11 dispCV: Coefficient of variation for dispersal
#  12 WhenBelow: When below local threshold (count in total or don't count or assume dead)
#  13 corrwith: Within population correlation - 
#    1 - (F, S, K Correlated)
#    2 - (F and S correlated, K independent)
#    3 - (F, S, K independent
#  14 DispersalDependsOnTargetPopK: Just like the variable name implie; boolean
#  15 DDBasis: Density basis type (char)
#  16 PopSpecificDD: Density dependence type is population specific (boolean - Yes/No)
#  17 DDforAllPop: Density dependence type for all populations (char)
#  18 UserDllFileName: Filename for user-defined density dependence function
#  19 TimeStepSize: Time step size
#  20 TimeStepUnits: Time step units 
#  21 SexStructure: Description of Sex Structure (char)
#  22 FemaleStages: Number of female stages
#  23 MatingSystem: Description of Mating system selected (char)
#  24 FemalesPerMale: Number of females per male as given by user
#  25 MalesPerFemale: Number of males per female as given by user
#  26 CVerror: sampling error for N
#  27 FirstTforRisk: Initial number of time steps to exclude from risk calculations
#  Population Parameters
#  28 PopList: Individual population data - saved as a 'nested list' structure, one 'list' for
#    each population.
#   Dispersal (Migration) Parameters
#   29 UseDispDistFunc: True if dispersal rates are based on dispersal distance function; false if
#     they are specified in the dispersal matrix
#   30 DispDistFunc: Dispersal-distance function parameters - a, b, c, Dmax - Mij = a exp(-Dij^c/b)
#  31 DispMatr: Dispersal matrix, as defined or calculated in 'fill.matrix.df.r'
#  Correlation Parameters
#   32 UseCorrDistFunc: True if correlations between populations is based on correlation distance
#     function; False if they are specified in the correlation matrix
#   33 CorrDistFunc: Correlation-distance function parameters - a, b, c - Cij = a exp(-Dij^c/b)
#   34 CorrMatr: User specified correlation matrix, as defined or calculated in 'fill.matrix.df.r'
#  Stage Matrices Information
#   35 StMatrNumber: number of stage matrices defined by user in *.mp file
#  36 StMatr: a 'list' object containing information about each stage matrix.  For each stage
#    matrix the following information is stored:
#    StMatrName: Name of the stage matrix
#    StMatrSurvMult: Survival multiplier for stage matrix (for generating new matrices; not used in simulation)
#    StMatrFecMult: Fecundity multiplier for stage matrix (for generating new matrices; not used in simulation)
#    Matr: a matrix object with numeric values
#  Standard Deviation Matrices Information
#  37 SDMatrNumber: number of st. dev. matrices defined by user in *.mp file
#  38 SDMatr: a 'list' object containing information for each st.dev. matrix.  For each st.dev.
#    matrix the following information is stored:
#    SDMatrName: Name of the st.dev. matrix
#    Matr: a matrix object of numeric values
#  39 ConstraintsMatr: Constraints Matrix
#   40 StMig: Stage Relative Migration (Dispersal) Rates
#   41 Cat1EffMat: Catastrophe 1 effects on vital rates; one row per stage, seperated by spaces
#   42 Cat1EffNst: Catastrophe 1 effects on abundances
#   43 Cat2EffMat: Catastrophe 2 effects on vital rates; one row per stage, seperated by spaces
#   44 Cat2EffNst: Catastrophe 2 effects on abundances
#   45 StInit: Initial Abundance for each stage for each population
#   Stage Properties
#   46 StProp: a 'list' object containing the following information for each stage
#     StName: Name of stage
#     StWeight: Relative weight of stage
#    StExclude: Exclude from total (boolean: TRUE or FALSE)
#    StBasisForDD: Basis for density-dependence (boolean: TRUE or FALSE)
#    StBreeding: Proportion of stage physically capable of breeding
#   Population Management Actions
#  47 NPopManage: Number of population managment actions
#   48 PopManageProp: All of the lines associated with population management, unpartitioned
#  49 ExtinctThr:  Extinction threshold
#  50 ExplodeThr:  Explosion threshold
#  51 stepsize: step size
#  52 PopData_df: Population level information in data.frame format
###################################################################################################

# Inform the user that the mp.read() function has been called.
print( paste( "Begin mp.read function with file: ", mpFile ) )
# Save the file path to mpFile
mpFilePath <- mpFile
# Read *.mp files into a long unsorted 'list' structure, 1 element per line in the file.
mpFile <- readLines(mpFile)
# Clear first line of "map=\xff". This step is needed because the "\x" leads to regex problems :(
mpFile[1] <- sub('map.*','',mpFile[1])

# Set values that will be used through out the function
FirstPopLine <- 45 # Line 45 correlates with first population data line
MigrationLine <- which(mpFile == "Migration")
CorrLine <- which(mpFile == "Correlation")
ConstraintsLine <- which(mpFile == "Constraints Matrix")

# Get Metapop Version information using metapopversion() function
metaVer <- metapopversion(mpFile)

# Create mp.file list and initiate with length 0
mp.file <- vector("list",length = 0)

# MaxRep: Number of replications (integer; 1-100000)
mp.file$MaxRep <- as.numeric(mpFile[7])

# MaxDur: Duration of simulation (integer; 1-500)
mp.file$MaxDur <- as.numeric(mpFile[8])

# Demog_Stoch: Use demographic stochasticity? (boolean)
mp.file$Demog_Stoch <- as.logical(mpFile[9])

# Stages: Number of stages (integer; 1-50)
# Stages.IgnoreConstraints: Ingore constraints? (boolean - TRUE/FALSE)
stageLine <- unlist(strsplit(mpFile[10],' '))
mp.file$Stages <- as.numeric(stageLine[1])
mp.file$Stages.IgnoreConstraints <- as.logical(stageLine[2])

# Cat1: Information associated with catastrophe 1
#  This information is subdivided in RAMAS Input, but will not be subdivided
#  in the first version of the sensitivity analysis package
mp.file$Cat1 <- mpFile[11:17]

# Cat2: Information associated with catastrophe 2
mp.file$Cat2 <- mpFile[18:25]

# DDActing: Information regarding density dependence
mp.file$DDActing <- mpFile[26]

# Distrib: Distribution (Normal or Lognormal) to use for Environmental
#  Stochasticity.
# AdvancedStochSetting: Advanced Environmental Stoch. Settings: (0 or PoolVar or 
#  NegSurvCorr)
stochLine <- unlist(strsplit(mpFile[27],','))
mp.file$Distrib <- stochLine[1]
mp.file$AdvancedStochSetting <- stochLine[2]

# dispCV: Coefficient of variation for dispersal
mp.file$dispCV <- as.numeric(mpFile[28])

# WhenBelow: When below local threshold (count in total or don't count or assume dead)
mp.file$WhenBelow <- mpFile[29]

# corrwith: Within population correlation - 
#  1 - (F, S, K Correlated)
#  2 - (F and S correlated, K independent)
#  3 - (F, S, K independent
mp.file$corrwith <- mpFile[30]

# DispersalDependsOnTargetPopK: Just like the variable name implie; boolean
mp.file$DispersalDependsOnTargetPopK <- mpFile[31]

# DDBasis: Density basis type (char). Possiblities are 'AllStages', 'SelectedStages', 'FecundityWeighted'
mp.file$DDBasis <- mpFile[32]

# PopSpecificDD: Density dependence type is population specific (boolean - Yes/No)
mp.file$PopSpecificDD <- mpFile[33]

# DDforAllPop: Density dependence type for all populations (char)
mp.file$DDforAllPop <- mpFile[34]

# UserDllFileName: Filename for user-defined density dependence function
mp.file$UserDllFileName <- mpFile[35]

# TimeStepSize: Time step size (integer)
mp.file$TimeStepSize <- mpFile[36]

# TimeStepUnits: Time step units (char)
mp.file$TimeStepUnits <- mpFile[37]

# SexStructure: Description of Sex Structure (char)
mp.file$SexStructure <- mpFile[38]

# FemaleStages: Number of female stages (integer)
mp.file$FemaleStages <- mpFile[39]

# MatingSystem: Description of Mating system selected (char)
mp.file$MatingSystem <- mpFile[40]

# FemalesPerMale: Number of females per male as given by user (number)
mp.file$FemalesPerMale <- as.numeric(mpFile[41])

# MalesPerFemale: Number of males per female as given by user
mp.file$MalesPerFemale <- as.numeric(mpFile[42])

# CVerror: sampling error for N
mp.file$CVerror <- as.numeric(mpFile[43])

# FirstTforRisk: Initial number of time steps to exclude from risk calculations
mp.file$FirstTforRisk <- as.numeric(mpFile[44])

# ----------------------------------------------------------------------------------------------- #
# PopList: Population level information 
print( "mp.read: Reading population information")
# First determine the number of populations (popNumber).  In version 5.0, population data begins at line 45
# and the 'Migration' section begins immediately after the last population's data
# The population level information is stored in two different structures - 1) a List format that 
# corresponds with how information is stored in RAMAS and 2) a data.frame that is convenient
# for functions.
##### TO DO - Below is setup for if there is no user defined variables
PopNumber <- MigrationLine - FirstPopLine
if (PopNumber < 1) {
  stop("Error: mp.read: Insufficient number of populations. Check 'first pop. line' value in mp.read.r")
}
PopRawData <- mpFile[FirstPopLine:(FirstPopLine + PopNumber - 1)]
# Create population data list object
PopData <- vector('list',length=0) # Initiate a PopData list
AllPopData <- vector('list',length=0) # Initiate a AllPopData list
PopData_df_rownames <- vector() # Initiate a row names vector
if ( metaVer == 50 ) {
  PopData_df_rownames <- c("name","X_coord","Y_coord","InitAbund","DensDep","MaxR","K","Ksdstr","Allee","KchangeSt","DD_Migr","Cat1.Multiplier","Cat1.Prob","IncludeInSum","StageMatType","RelFec","RelSur","localthr","Cat2.Multiplier","Cat2.Prob","SDMatType","TargetPopK","Cat1.TimeSinceLast","Cat2.TimeSinceLast","RelDisp")
  for ( pop in 1:PopNumber ) {
    popLine <- unlist(strsplit( PopRawData[ pop ], ',' ))
    PopData$name <- popLine[1]  
    PopData$X_coord <- as.numeric(popLine[2])
    PopData$Y_coord <- as.numeric(popLine[3])
    PopData$InitAbund <- popLine[4]
    PopData$DensDep <- popLine[5]
    PopData$MaxR <- popLine[6]
    PopData$K <- popLine[7]
    PopData$Ksdstr <- popLine[8]
    PopData$Allee <- popLine[9]
    PopData$KchangeSt <- popLine[10]
    PopData$DD_Migr <- popLine[11]
    PopData$Cat1.Multiplier <- popLine[12]
    PopData$Cat1.Prob <- popLine[13]
    PopData$IncludeInSum <- popLine[14]
    PopData$StageMatType <- popLine[15]
    PopData$RelFec <- popLine[16]
    PopData$RelSur <- popLine[17]
    PopData$localthr <- popLine[18]
    PopData$Cat2.Multiplier <- popLine[19]
    PopData$Cat2.Prob <- popLine[20]
    PopData$SDMatType <- popLine[21]
    PopData$TargetPopK <- popLine[22]
    PopData$Cat1.TimeSinceLast <- popLine[23]
    PopData$Cat2.TimeSinceLast <- popLine[24]
    PopData$RelDisp <- popLine[25]
    AllPopData[[pop]] <- PopData # Add PopData to full list of population data
  }
} else if ( metaVer >= 51 ) {
  PopData_df_rownames <- c("name","X_coord","Y_coord","InitAbund","DensDep","MaxR","K","Ksdstr","Allee","KchangeSt","DD_Migr","Cat1.Multiplier","Cat1.Prob","IncludeInSum","StageMatType","RelFec","RelSur","localthr","Cat2.Multiplier","Cat2.Prob","SDMatType","TargetPopK","Cat1.TimeSinceLast","Cat2.TimeSinceLast","RelDisp","RelVarFec","RelVarSurv")
  for ( pop in 1:PopNumber ) {
    popLine <- unlist(strsplit( PopRawData[ pop ], ',' ))
    PopData$name <- popLine[1]  
    PopData$X_coord <- as.numeric(popLine[2])
    PopData$Y_coord <- as.numeric(popLine[3])
    PopData$InitAbund <- as.numeric(popLine[4])
    PopData$DensDep <- popLine[5]
    PopData$MaxR <- as.numeric(popLine[6])
    PopData$K <- as.numeric(popLine[7])
    PopData$Ksdstr <- as.numeric(popLine[8])
    PopData$Allee <- as.numeric(popLine[9])
    PopData$KchangeSt <- popLine[10]
    PopData$DD_Migr <- popLine[11]
    PopData$Cat1.Multiplier <- popLine[12]
    PopData$Cat1.Prob <- popLine[13]
    PopData$IncludeInSum <- popLine[14]
    PopData$StageMatType <- popLine[15]
    PopData$RelFec <- popLine[16]
    PopData$RelSur <- popLine[17]
    PopData$localthr <- popLine[18]
    PopData$Cat2.Multiplier <- popLine[19]
    PopData$Cat2.Prob <- popLine[20]
    PopData$SDMatType <- popLine[21]
    PopData$TargetPopK <- popLine[22]
    PopData$Cat1.TimeSinceLast <- popLine[23]
    PopData$Cat2.TimeSinceLast <- popLine[24]
    PopData$RelDisp <- popLine[25]
    PopData$RelVarFec <- popLine[26]
    PopData$RelVarSurv <- popLine[27]
#     ## Look for user defined density dependence parameters
#     if (length(popLine)>27){
#       for (udd_par in 1:(length(popLine)-27)){
#         PopData$UserDef
#       }
#     }
    AllPopData[[pop]] <- PopData # Add PopData to full list of population data
  }
}
mp.file$PopList <- AllPopData
# Create population data data.frame
PopData_df <- read.csv( mpFilePath, header=FALSE, skip=44, nrows=PopNumber )
###browser()
### Only use as many columns as there are elements in the PopData_df_rownames.  The rest of the columns
### are associated with user defined density dependence parameters, not used at this time.
### PopData_df <- PopData_df[1:length(PopData_df_rownames)]
# Turn NAs into empty strings
PopData_df[is.na(PopData_df)] <- ''
# Check if PopData_df includes user defined d-d values
if (ncol(PopData_df)>27){
  # Get the number of userd defined d-d pars
  Num_udd_pars <- ncol(PopData_df)-27
  udd_names <- paste('udd_',1:Num_udd_pars,sep='')
  PopData_df_rownames <- c(PopData_df_rownames,udd_names)
}
# Assign columns names
names(PopData_df) <- PopData_df_rownames
###browser()
# ----------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------- #
# Dispersal (Migration) Data 
print( "mp.read: Reading dispersal (migration) information" )
# UseDispDistFunc: True if dispersal rates are based on dispersal distance function; false if
#   they are specified in the dispersal matrix
mp.file$UseDispDistFunc <- as.logical( mpFile[MigrationLine + 1] )
# DispDistFunc: Dispersal-distance function parameters - a, b, c, Dmax - Mij = a exp(-Dij^c/b)
mp.file$DispDistFunc <- as.numeric(unlist(strsplit( mpFile[MigrationLine + 2],',' )))
# DispMatr: User specified dispersal matrix. If the user selected dispersal based on 
#   disp-dist function then there are no rows for dispersal matrix. The definition of 
#  DispMatr has non-intuitive indexing to account for the fact that this first two lines
#  after "Migration" are not part of the matrix, they are UseDispDistFunc and
#  DispDistFunc, respectively.
if ( mp.file$UseDispDistFunc ) {
  # migLines is a variable used to identify the number of lines used to define the
  # disp-dist func parameters and the disp matrix, if one is defined
  migLines <- 1 #Only one line necessary for migration parameters
  mp.file$DispMatr <- fill.matrix.df( PopData_df, mp.file$DispDistFunc, 'disp' )
} else {
  migLines <- (1 + PopNumber)
  DispMatr <- mpFile[ (MigrationLine + 3):(MigrationLine + 1 + migLines) ]
  mp.file$DispMatr <- matrix(as.numeric(unlist(strsplit( DispMatr,',' ))), nrow = PopNumber, byrow = TRUE)
}
# ----------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------- #
# Correlation Data
print( "mp.read: Reading correlation information")
# UseCorrDistFunc: True if correlations between populations is based on correlation distance
#   function; False if they are specified in the correlation matrix
mp.file$UseCorrDistFunc <- as.logical( mpFile[CorrLine +1] )
# CorrDistFunc: Correlation-distance function parameters - a, b, c - Cij = a exp(-Dij^c/b)
mp.file$CorrDistFunc <- as.numeric(unlist(strsplit( mpFile[CorrLine + 2],',' )))
# CorrMatr: User specified correlation matrix. Similar to dispersl matrix (See notes for DispMatr 
# above), however the matrix is a lower triangular matrix, thus it is not saved as a matrix object
if ( mp.file$UseCorrDistFunc ) {
  corrLines <- 1
  # Create a Correlation Matrix from the Correlation Distance function.  If there is only one 
  # population, then the number 1 is returned
  mp.file$CorrMatr <- fill.matrix.df( PopData_df, mp.file$CorrDistFunc, 'corr' )
} else {
  print("Using Correlation matrix") ### WARNING LINE
  # Define a new function used to read correlation distance matrices
  # addZeroes: used to read a correlation distance matrix.  This function adds zeroes to each line
  # of the lower triangular corr-dist matrix stored in the *.mp file
  addZeroes <- function( vect ) { 
    zeroes <- PopNumber - length(vect)
    new_vect <- c( vect, rep(0, zeroes) )
    return(new_vect)
  }
  corrLines <- (1 + PopNumber)
  # If user defined, the corr-matr is a lower-triangular matrix. To make comparisons easy, the 
  # matrix is transformed into a 'matrix' object
  #
  # Read corrMatr as list from readLines input
  corrMatr <- mpFile[ (CorrLine+3):(CorrLine+2+PopNumber) ] 
  # Split each list element by comma
  corrMatr <- strsplit( corrMatr, ',' )
  # Make each list element a vector of numeric values
  corrMatr <- lapply( corrMatr, as.numeric)
  # Apply the addZeroes function
  corrMatr <- lapply( corrMatr, addZeroes )
  # Make the 'list' into a 'matrix'
  corrMatr <- do.call( rbind, corrMatr )
  # Set CorrMatr in mp.file 'list'
  mp.file$CorrMatr <- corrMatr
}
# ----------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------- #
# Stage Matrices Information
# Determine the number of stage matrix types (i.e., how many stage matrices are defined)
# The number of stage matices (stage matrix types) is the first number of the line that 
# contains the phrase "stage matrix" in it.
stgMatrLine <- grep('stage matrix', mpFile, ignore.case=TRUE)
StMatrNumber <- unlist(strsplit( mpFile[stgMatrLine],' ' ))
StMatrNumber <- as.numeric( StMatrNumber[1] )
# StMatrNumber: number of stage matrices defined by user in *.mp file
mp.file$StMatrNumber <- StMatrNumber
#
# StMatr: a list object containing information about each stage matrix.
StMatr <- vector('list',length=0) # Create an empty stage matrix
oneMatr <- 4 + mp.file$Stages # The number of lines of information for one matrix
allMatr <- StMatrNumber * oneMatr # The number of lines of information for all matrices
# Extract from mpFile all of the lines of information for all of the matrices
AllStMatrLines <- mpFile[ (stgMatrLine + 1):(stgMatrLine + 1 + allMatr) ]
#
AllStMatr <- vector('list',length=0) # Create an empty list of stage matrices
for ( matr in 1:mp.file$StMatrNumber ) {
  LineAdd <- matr - 1
  # StMatrName: Name of the stage matrix
  StMatr$StMatrName <- AllStMatrLines[ 1 + (LineAdd*oneMatr) ]
  # StMatrSurvMult: Survival multiplier for stage matrix (for generating new matrices; not used in simulation)
  StMatr$StMatrSurvMult <- AllStMatrLines[ 2 + (LineAdd*oneMatr) ]
  # StMatrFecMult: Fecundity multiplier for stage matrix (for generating new matrices; not used in simulation)
  StMatr$StMatrFecMult <- AllStMatrLines[ 3 + (LineAdd*oneMatr) ]
  # Matr: a matrix object with numeric values
  dumbMatr <- AllStMatrLines[ (5 + (LineAdd*oneMatr)):( (4 + mp.file$Stages) + (LineAdd*oneMatr) ) ]
  StMatr$Matr <- matrix(as.numeric(unlist(strsplit( dumbMatr,' ' ))), nrow = mp.file$Stages, byrow = TRUE)
  AllStMatr[[matr]] <- StMatr
}
mp.file$StMatr <- AllStMatr
# ----------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------- #
# Standard Deviation Matrices Information
# Determine the number of std. dev. matrix types.
# The number of standard deviation matrices (st dev matrix types) is the first number of the line that 
# contains the phrase "st.dev. matrix" in it.
sdMatrLine <- grep('st.dev. matrix', mpFile, ignore.case=TRUE)
sdMatrNumber <- unlist(strsplit( mpFile[sdMatrLine],' ' ))
sdMatrNumber <- as.numeric( sdMatrNumber[1] )
# SDMatrNumber: number of st. dev. matrices defined by user in *.mp file
mp.file$SDMatrNumber <- sdMatrNumber
#
# SDMatr: a list object containing information for each st.dev. matrix.
SDMatr <- vector('list',length=0) # Create an empty st.dev. matrix list object
oneSDMatr <- 1 + mp.file$Stages # The number of lines of information for one st.dev. matrix
allSDMatr <- sdMatrNumber * oneSDMatr # The number of lines of information for all st.dev. matrices
# Extract from mpFile all of the lines of information for all of the matrices
AllSDMatrLines <- mpFile[ (sdMatrLine + 1):(sdMatrLine + 1 + allSDMatr) ]
#
AllSDMatr <- vector('list',length=0) # Create an empty list of std deviation matrices
for ( matr in 1:mp.file$SDMatrNumber ) {
  LineAdd <- matr - 1 # Addition factor to skip to the lines associated with the matrix of interest
  # SDMatrName: Name of the st.dev. matrix
  SDMatr$SDMatrName <- AllSDMatrLines[ 1 + (LineAdd*oneSDMatr) ]
  # Matr: a matrix object of numeric values
  dumbMatr <- AllSDMatrLines[ (2 + (LineAdd*oneSDMatr)):( (1 + mp.file$Stages) + (LineAdd*oneSDMatr) ) ]
  SDMatr$Matr <- matrix(as.numeric(unlist(strsplit( dumbMatr,' ' ))), nrow = mp.file$Stages, byrow = TRUE)
  AllSDMatr[[matr]] <- SDMatr
}
mp.file$SDMatr <- AllSDMatr
# ----------------------------------------------------------------------------------------------- #

# ConstraintsMatr: Constraints Matrix
dumbMatr <- mpFile[ (ConstraintsLine + 1):(ConstraintsLine + mp.file$Stages) ]
mp.file$ConstraintsMatr <- matrix(as.numeric(unlist(strsplit( dumbMatr,' '))),nrow = mp.file$Stages, byrow = TRUE)

# StMig: Stage Relative Migration (Dispersal) Rates
mp.file$StMig <- as.numeric(unlist(strsplit( mpFile[ ConstraintsLine + mp.file$Stages + 1 ],' ' )))

# Cat1EffMat: Catastrophe 1 effects on vital rates; one row per stage, seperated by spaces
dumbMatr <- mpFile[ (ConstraintsLine + mp.file$Stages + 2):(ConstraintsLine + 2*mp.file$Stages + 1) ]
mp.file$Cat1EffMat <- matrix(as.numeric(unlist(strsplit( dumbMatr,' ' ))), nrow = mp.file$Stages, byrow = TRUE)

# Cat1EffNst: Catastrophe 1 effects on abundances
mp.file$Cat1EffNst <- as.numeric(unlist(strsplit( mpFile[ ConstraintsLine + 2*mp.file$Stages + 2 ],' ' )))

# Cat2EffMat: Catastrophe 2 effects on vital rates; one row per stage, seperated by spaces
dumbMatr <- mpFile[ (ConstraintsLine + 2*mp.file$Stages + 3):(ConstraintsLine + 3*mp.file$Stages + 2) ]
mp.file$Cat2EffMat <- matrix(as.numeric(unlist(strsplit( dumbMatr,' ' ))), nrow = mp.file$Stages, byrow = TRUE)

# Cat2EffNst: Catastrophe 2 effects on abundances
mp.file$Cat2EffNst <- as.numeric(unlist(strsplit( mpFile[ ConstraintsLine + 3*mp.file$Stages + 3 ],' ' )))

# StInit: Initial Abundance for each stage for each population
InitAbLine <- ConstraintsLine + 3*mp.file$Stages + 4 # First line
InitAb <- mpFile[ InitAbLine:(InitAbLine + PopNumber - 1) ]
mp.file$StInit <- matrix(as.numeric(unlist(strsplit( InitAb,' ' ))), nrow = PopNumber, byrow = TRUE)

# ----------------------------------------------------------------------------------------------- #
# Stage Properties
InitStPropLine <- InitAbLine + PopNumber # First line of stage properties in the *.mp file
# StProp: a list object containing information for each stage
StProp <- vector('list',length=0) # Create an empty stage property list object
# For each stage there are five properties, thus five lines of information
allStProp <- 5 * mp.file$Stages
AllStPropLines <- mpFile[ InitStPropLine:(InitStPropLine + allStProp - 1) ]
#
AllStProp <- vector('list',length=0) # Create an empty list of St.Prop. lists
for ( stg in 1:mp.file$Stages ) {
  firstinfo <- (stg - 1)*5 + 1 # First line of info for stage 'stg'
  # StName: Name of stage
  StProp$StName <- AllStPropLines[ firstinfo ]
  # StWeight: Relative weight of stage
  StProp$StWeight <- AllStPropLines[ firstinfo + 1 ]
  #StExclude: Exclude from total (boolean: TRUE or FALSE)
  StProp$StExclude <- AllStPropLines[ firstinfo + 2 ]
  #StBasisForDD: Basis for density-dependence (boolean: TRUE or FALSE)
  StProp$StBasisForDD <- AllStPropLines[ firstinfo + 3 ]
  #StBreeding: ###### WHAT IS THIS????
  StProp$StBreeding <- AllStPropLines[ firstinfo + 4 ]
  AllStProp[[ stg ]] <- StProp
}
mp.file$StProp <- AllStProp
# ----------------------------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------------------------- #
# Population Management Actions
# NPopManage: The number of population managment actions
mgmntLine <- grep('pop mgmnt',mpFile,ignore.case=TRUE)
NPopManage <- unlist(strsplit(mpFile[mgmntLine],' '))
NPopManage <- as.numeric(NPopManage[1])
mp.file$NPopManage <- NPopManage
if ( NPopManage > 0 ) {
  # PopManageProp: All of the lines associated with population management, unpartitioned
  ####mp.file$PopManageProp <- mpFile[ (mgmntLine + 1):(mgmntLine + NPopManage) ]
  PopManageProp.colNames <- c('active', 'mng.type', 'from.pop', 'to.pop', 
                              'begin.time','end.time','period','when','num.or.prop',
                              'number','proportion',
                              'from.stage','to.stage','cond.type','cond.abund.low',
                              'cond.abund.high','from.all.stgs','cond.quant.2','cond.func.N1',
                              'cond.func.N2','abund.each.stg_div_all.stg',
                              'abund.each.pop_div_all.pop')
    
  mp.file$PopManageProp <- read.table(  mpFilePath , skip=mgmntLine, nrows=NPopManage,
                                        col.names=PopManageProp.colNames)
} else {
  mp.file$PopManageProp <- "NA" # Fill with 'NA' value
}
# ----------------------------------------------------------------------------------------------- #
# The next three elements of the list are defined based on the position of the population
# management line. This may change in the future and have to be adjusted

# ExtinctThr: Extinction Threshold
mp.file$ExtinctThr <- as.numeric( mpFile[ mgmntLine + NPopManage + 1 ] )

# ExplodeThr: Explosion Threshold
mp.file$ExplodeThr <- as.numeric( mpFile[ mgmntLine + NPopManage + 2 ] )

# stepsize: Stepsize
mp.file$stepsize <- as.numeric( mpFile[ mgmntLine + NPopManage + 3 ] )

# ----------------------------------------------------------------------------------------------- #
# The last element of the list is the population data data.frame
mp.file$PopData_df <- PopData_df

return( list( version = metaVer, mp.file = mp.file) )
} # End mp.read function

