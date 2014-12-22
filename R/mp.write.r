mp.write <- function( mp.new, version, mp.new.file ) {
# Purpose: Write a new *.mp file. 
# 
# Args: 
#  mp.new: mp.file list structure, the result of the mp.read.r or sensitivity.r functions
#
#  version: the Metapop version number of the mp.file
#
#  mp.new.file: file name for the new *.mp file to be written
#
# Returns:
#  None within R
#
# The input mp.new structure should be the structure resulting from either the mp.read.r or the
# sensitivity.r scripts.  This script produces a file readable by the RAMAS Metapop program.
###################################################################################################

# Inform the user that the mp.write() function has been called.
print( paste( "Begin mp.write function with file: ", mp.new.file ) )
# Remove any files with the same name as the new *.mp file
unlink(mp.new.file)

# First line - includes version number and skips five lines
version.out <- paste( substr(version, 1, 1), ".", substr(version, 2, 2), sep = "")
line.1 <- paste("Metapopulation input file (", version.out, ") map=", "\n\n\n\n\n", sep = "")
write( line.1, file = mp.new.file, append = TRUE )

# Number of replications
write( mp.new$MaxRep, file = mp.new.file, append = TRUE )

# Duration of simulation
write( mp.new$MaxDur, file = mp.new.file, append = TRUE )

# Use demographic stochasticity?
write( mp.new$Demog_Stoch, file = mp.new.file, append = TRUE )

# Number of stages and Ignore constraints?
line.10 <- paste( mp.new$Stages, mp.new$Stages.IgnoreConstraints, sep = " ")
write( line.10, file = mp.new.file, append = TRUE )

# Catastrophe-1 information 
write( mp.new$Cat1, file = mp.new.file, append = TRUE )
# Catastrophe-2 information 
write( mp.new$Cat2, file = mp.new.file, append = TRUE )

# Density Dependence Acting
write( mp.new$DDActing, file = mp.new.file, append = TRUE )

# Enviromental stochasticity distribution and advanced stochasticity settings
line.27 <- paste( mp.new$Distrib, mp.new$AdvancedStochSetting, sep = "," )
write( line.27, file = mp.new.file, append = TRUE )

# Coefficient of Variation for Dispersal
write( mp.new$dispCV, file = mp.new.file, append = TRUE )

# When below local threshold?
write( mp.new$WhenBelow, file = mp.new.file, append = TRUE )

# Within population correllation
write( mp.new$corrwith, file = mp.new.file, append = TRUE )

# Dispersal depends on target population K
write( mp.new$DispersalDependsOnTargetPopK, file = mp.new.file, append = TRUE )

# Density dependence basis type
write( mp.new$DDBasis, file = mp.new.file, append = TRUE )

# Density type is population specific?
write( mp.new$PopSpecificDD, file = mp.new.file, append = TRUE )

# Density dependence type for all populations
write( mp.new$DDforAllPop, file = mp.new.file, append = TRUE )

# File name for user defined density dependence function
write( mp.new$UserDllFileName, file = mp.new.file, append = TRUE )

# Time step size
write( mp.new$TimeStepSize, file = mp.new.file, append = TRUE )

# Time step units
write( mp.new$TimeStepUnits, file = mp.new.file, append = TRUE )

# Sex structure
write( mp.new$SexStructure, file = mp.new.file, append = TRUE )

# Number of female stages
write( mp.new$FemaleStages, file = mp.new.file, append = TRUE )

# Mating system
write( mp.new$MatingSystem, file = mp.new.file, append = TRUE )

# Females per male
write( mp.new$FemalesPerMale, file = mp.new.file, append = TRUE )

# Males per female
write( mp.new$MalesPerFemale, file = mp.new.file, append = TRUE )

# Sampling error for N
write( mp.new$CVerror, file = mp.new.file, append = TRUE )

# Initial number of time steps to exclude from risk calculations
write( mp.new$FirstTforRisk, file = mp.new.file, append = TRUE )

# Population Specific parameters
print( "mp.write: Writing population specific parameters")
# If the pop. parameters are not in a data.frame format, take them from mp.file$PopData_df
if ( !is.data.frame( mp.new$PopList ) ) {
  print("WARNING: Writing populaiton specific parameters from PopData_df.  If you made changes to PopList elements, and did not also change the PopData_df elements, your changes may not be written to a new *.mp file.")
  mp.new$PopList <- mp.new$PopData_df 
} 
write.table( mp.new$PopList, file = mp.new.file, append=TRUE, sep = ",", row.names=FALSE, col.names=FALSE, quote=FALSE )
# Define a variable that is the number of populaitons
pop.num <- nrow(mp.new$PopList)

# Migration parameters
print( "mp.write: Writing migration parameters" )
write( "Migration", file = mp.new.file, append = TRUE )
# Use dispersal distance function?
write( mp.new$UseDispDistFunc, file = mp.new.file, append = TRUE )
# Dispersal distance function parameters
write( mp.new$DispDistFunc, file = mp.new.file, append = TRUE, sep = ",")
# If UseDispDistFunc = FALSE then write the dispersal matrix in full
if ( !mp.new$UseDispDistFunc ) {
  # Dispersal matrix
  write.table( mp.new$DispMatr, file = mp.new.file, append=TRUE, sep = ",", eol = ",\n", row.names=FALSE, col.names=FALSE, quote=FALSE )
}

# Correlation paramters
print( "mp.write: Writing correlation parameters" )
write( "Correlation", file = mp.new.file, append = TRUE )
# Use correlation distance function?
write( mp.new$UseCorrDistFunc, file = mp.new.file, append = TRUE )
# Correlation distance function parameters
write( mp.new$CorrDistFunc, file = mp.new.file, append = TRUE, sep = "," )
# Correlation matrix
# If UseCorrDistFunc = FALSE then write the correlation matrix in full
if ( !mp.new$UseCorrDistFunc ) {
  # This writes a full matrix to the *.mp file, which is not the correct format.  However, the RAMAS
  # program appears to only consider the part of the matrix that is appropriate, and when the *mp file
  # is saved in RAMAS, it has the correct format.
  write.table( mp.new$CorrMatr, file = mp.new.file, append=TRUE, sep = ",", eol = ",\n", row.names=FALSE, col.names=FALSE, quote=FALSE )
}


####################### KTS: changed to accommodate multiple stage matrices: 8/30/2013

# Stage matrix parameters
templine <- paste(mp.new$StMatrNumber," type(s) of stage matrix", sep = "")   
write( templine, file = mp.new.file, append = TRUE )

for(matr in 1: mp.new$StMatrNumber){
	# Matrix name
	write( mp.new$StMatr[[matr]]$StMatrName, file = mp.new.file, append = TRUE )
	# Survival multiplier
	write( mp.new$StMatr[[matr]]$StMatrSurvMult, file = mp.new.file, append = TRUE )
	# Fecundity multiplier
	write( mp.new$StMatr[[matr]]$StMatrFecMult, file = mp.new.file, append = TRUE )
	# Zero line
	write( 0, file = mp.new.file, append = TRUE )
	# Stage matrix
	for ( stage in 1:mp.new$Stages ) {
	  write( mp.new$StMatr[[matr]]$Matr[stage, ], file = mp.new.file, append = TRUE, ncolumns = mp.new$Stages )
	}
}

# Standard deviation matrix parameters
templine <- paste(mp.new$SDMatrNumber," type(s) of st.dev. matrix", sep = "")   
write( templine, file = mp.new.file, append = TRUE )

for(matr in 1:mp.new$SDMatrNumber){
	# SD Matrix name
	write( mp.new$SDMatr[[matr]]$SDMatrName, file = mp.new.file, append = TRUE )
	# SD Matrix
	for ( stage in 1:mp.new$Stages ) {
	  write( mp.new$SDMatr[[matr]]$Matr[stage, ], file = mp.new.file, append = TRUE, ncolumns = mp.new$Stages )
	}
}

########################### END KTS CHANGES

# Constraints matrix
write( "Constraints Matrix", file = mp.new.file, append = TRUE )
for ( stage in 1:mp.new$Stages ) {
  write( mp.new$ConstraintsMatr[stage, ], file = mp.new.file, append = TRUE, ncolumns = mp.new$Stages )
}

# StMig: Stage Relative Migration (Dispersal) Rates
write( mp.new$StMig, file = mp.new.file, append = TRUE, ncolumns = mp.new$Stages )

# Catastrophe effects
# Cat1EffMat: Catastrophe 1 effects on vital rates; one row per stage, seperated by spaces
for ( stage in 1:mp.new$Stages ) {
  write( mp.new$Cat1EffMat[stage, ], file = mp.new.file, append = TRUE, ncolumns = mp.new$Stages )
}
# Cat1EffNst: Catastrophe 1 effects on abundances
write( mp.new$Cat1EffNst, file = mp.new.file, append = TRUE, ncolumns = mp.new$Stages )
# Cat2EffMat: Catastrophe 2 effects on vital rates; one row per stage, seperated by spaces
for ( stage in 1:mp.new$Stages ) {
  write( mp.new$Cat2EffMat[stage, ], file = mp.new.file, append = TRUE, ncolumns = mp.new$Stages )
}
# Cat2EffNst: Catastrophe 2 effects on abundances
write( mp.new$Cat2EffNst, file = mp.new.file, append = TRUE, ncolumns = mp.new$Stages )

# StInit: Initial Abundance for each stage for each population
for ( p in 1:pop.num ) {
  write( mp.new$StInit[p, ], file = mp.new.file, append = TRUE, ncolumns = mp.new$Stages )
}

# Stage Properties
for ( stage in 1:mp.new$Stages ) {
  # There are five properties for each stage (as of version 5)
  write( mp.new$StProp[[ stage ]]$StName, file = mp.new.file, append = TRUE ) # Stage name
  write( mp.new$StProp[[ stage ]]$StWeight, file = mp.new.file, append = TRUE ) # Relative weight
  write( mp.new$StProp[[ stage ]]$StExclude, file = mp.new.file, append = TRUE ) # Exclude from total?
  write( mp.new$StProp[[ stage ]]$StBasisForDD, file = mp.new.file, append = TRUE ) # Basis for Dens Dep
  write( mp.new$StProp[[ stage ]]$StBreeding, file = mp.new.file, append = TRUE ) # Breeding
}

# Number of population managment actions
NPopManage <- paste( mp.new$NPopManage, " (pop mgmnt)" )
write( NPopManage, file = mp.new.file, append = TRUE )

# Population management properties
if ( all( mp.new$PopManageProp != "NA" ) ) {
  write.table( mp.new$PopManageProp, file = mp.new.file, append=TRUE, sep = "\t", row.names=FALSE, col.names=FALSE, quote=FALSE )
}

# Extinction threshold
write( mp.new$ExtinctThr, file = mp.new.file, append = TRUE )

# Explosion threshold
write( mp.new$ExplodeThr, file = mp.new.file, append = TRUE )

# stepsize
write( mp.new$stepsize, file = mp.new.file, append = TRUE )

# End of file marker
write( "-End of file-", file = mp.new.file, append = TRUE )

} # End mp.write function
