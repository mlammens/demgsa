mp.results <- function( mpFile, mpList, spatial=FALSE, mac=FALSE, ptc=FALSE, ptcFiles='no file', 
  ptcFileIter='', habdyn=FALSE, hdhFile='no file' ) {
# Calculate various endpoints from the results of an MP simulation as carried out in RAMAS
# Use for Metapop version 5 and 5.1 file formats
# Based on ResOut.exe program written by R. Akcakaya
#
# Author: Matthew Aiello-Lammens
# Created: 9 January 2012
# Update: 21 February 2012; 5 March 2012 (added spatial section)
#
# Args:
#  mpFile: the full path name to the *.mp file to be examined
#
#  spatial: TRUE/FALSE - should spatial parameters be extracted from the *.mp file,
#    including factors such as dispersal, correlation, and factors included in
#    *.ptc files, if included
#
#  mac: TRUE/FALSE - Is the run on a mac or *nix machine vs a Windows machine
#
#  ptc: TRUE/FALSE - Is a *.ptc file included?  If yes, some spatial information can
#    be extracted from this file
#
#  ptcFiles: A single path to a ptc file or a vector of ptc files with which to extract
#    data from.
#
#  ptcFileIter: For each ptcFile, the corresponding iteration as noted in the HabDyn 
#    history file.  If no history file is used, iter values can be some user defined
#    values, however an error will be generated if the length of ptcfiles and
#    ptcFileIter do not match.
#
#  habdyn: TRUE / FALSE - Was HabDyn module used? Is there a corresponding HabDyn
#    history file?
#
#  hdhFile: The full path name to the HabDyn history (.txt) file
#
# Returns:
#  res.summ.df: A data frame of the endpoints calculated
#
###################################################################################################
# Inform the user that mp.results has begun
#print(paste('Begin mp.results funtion with file: ', mpFile))
  
### FRAL Specific Change - To analyze these data I need to have already read the mpFile

# # Check that *.mp file exists
# if ( !file.exists( mpFile ) ){
#   stop( paste( 'Did not find *.mp file:', mpFile) )
# }

# Read *.mp file using mp.read.results.r.  If the file does not contain results
# it will be reported in the read script.
#mp <- mp.read.results( mpFile )
mp <- mpList

# Get the *.mp parameter inputs
mp.in <- mp[1:52]
# Get the *.mp results
mp.res <- mp$results

# Created Results Summary (res.summ) list and initiate with length = 0
res.summ <- vector("list", length=0)

##### BEGIN REPORT OF INPUT VALUES (Indpendent Variables) #####
# ----------------------------------------------------------------------------------------------- #
# Metapopulation Initial Abundance: Use only values that are 'included in sum'
pop.initab <- mp.in$PopData_df$InitAbund
pop.includeinsum <- mp.in$PopData_df$IncludeInSum
res.summ$metapop.initab <- sum( pop.initab * pop.includeinsum )
# Calculate 50% decline threshold - to be used later in the script
Thr.50 <- res.summ$metapop.initab / 2
# ----------------------------------------------------------------------------------------------- #
# # Call GenTime.exe to determine various calculations of generation time and eigen values
# if( mac ) {
#   # WARNING: Running GenTime on Mac or Linux requires a Wine installation on the system
#   # and the call to the program can be temperamental depending on syntax used.
#   #wine <- '/Applications/Wine.app/Contents/Resources/bin/wine' # Works for my Mac
#   wine <- '/opt/local/bin/wine'
#   # Assuming GenTime.exe is in same directory as scripts
#   gen.exe <- paste( sens.base.dir, 'GenTime.exe ', sep='')
#   gen.call <- paste( wine, gen.exe , '"', mpFile , '"' )
#   gen <- system( gen.call, intern=TRUE, ignore.stderr = TRUE )
#   gen <- strsplit( gen, " ")
#   gen <- unlist( lapply( gen, as.numeric ) )
#   res.summ$Gen.Abar <- gen[1]
#   res.summ$Gen.Mu1 <- gen[2]
#   res.summ$Gen.Tgen <- gen[3]
#   res.summ$EigenVal <- gen[4]
# } else {
#   # Calling GenTime can be temperamental depending on syntax used.
#   gen.call <- paste( sens.base.dir,'GenTime ', '"', mpFile , '"',sep="" )
#   print( paste('GenTime.exe called as: ', gen.call) ) ### DEBUG LINE
#   gen <- system( gen.call, intern = TRUE )
#   gen <- strsplit( gen, " ")
#   gen <- unlist( lapply( gen, as.numeric ) )
#   res.summ$Gen.Abar <- gen[1]
#   res.summ$Gen.Mu1 <- gen[2]
#   res.summ$Gen.Tgen <- gen[3]
#   res.summ$EigenVal <- gen[4]
# }
## Use 'popbio' package functions to get generation time and eigen value
res.summ$EigenVal <- eigen.analysis(mp.in$StMatr[[1]]$Matr)$lambda1
res.summ$GenTime <- generation.time(mp$StMatr[[1]]$Matr)

# ----------------------------------------------------------------------------------------------- #
# Density Dependence Factors
#
# Determine denisty type:
# Is density dependence population specific? If yes, then determine dd.type from pop info,
# if no, then dd.type is DDforAllPop value
dd.pop.spec <- mp.in$PopSpecificDD
if ( dd.pop.spec == 'Yes' ) {
  dd.type.pop <- mp.in$PopData_df$DensDep
  # If there is more than one dd.type in this vector, then report dd.type as mixed
  if ( length(unique(dd.type.pop)) > 1 ) { 
    dd.type <- 'Mixed'
  } else { 
    dd.type <- unique(dd.type.pop)
  }
} else if ( dd.pop.spec == 'No' ) {
  dd.type <- mp.in$DDforAllPop
  # If dd.type is user defined, also include DLL file
  if ( dd.type == 'UD') {
    dd.type <- paste( dd.type, mp.in$UserDllFileName )
  }
}
res.summ$dd.type <- dd.type
# ----------------------------------------------------------------------------------------------- #
# Rmax values: Calculate the mean Rmax values over all poplations
rmax.pop <- mp.in$PopData_df$MaxR
res.summ$Rmax <- mean(rmax.pop)
# ----------------------------------------------------------------------------------------------- #
# Growth Rate: Based on density dependence type
#
# A vector of chars that define the different density types for which growth rate
# would equal the eigen value of the stage matrix
dd.eig <- c('EX','CE','EA','CA')
# Set growth rate according to dd.type recorded
if ( any( dd.type == dd.eig )) {
  res.summ$GrowthRt <- res.summ$EigenVal
} else {
  # Assume that if dd.type is not in dd.eig, then use Rmax.  This will include
  # user defined dd.type
  res.summ$GrowthRt <- res.summ$Rmax
}
# ----------------------------------------------------------------------------------------------- #
# Variability Metrics
#
# Get constraints matrix - 1 equals survival, 0 equals fecundity
cons <- mp.in$ConstraintsMatr
# Get st.dev matr
SD.matr <- mp.in$SDMatr[[1]]$Matr
# Get stage matrix
St.matr <- mp.in$StMatr[[1]]$Matr
# Seperate survival and fecundity elements of matrix
surv.matr <- St.matr * cons
fec.matr <- St.matr * (1-cons)
# Calculate a variance matrix
var.matr <- SD.matr^2
# Average fecundity st.dev of all non-zero F values
res.summ$fec.stdev.avg <- sqrt( mean( var.matr[ which( fec.matr > 0 ) ] ) )
# Average survival st.dev of all non-zero S values
res.summ$surv.stdev.avg <- sqrt( mean( var.matr[ which( surv.matr > 0 ) ] ) )
# Average standard deviation, averaged over all non-zero matrix elements
res.summ$stdev.avg <- sqrt( mean( var.matr[ which( St.matr > 0 ) ] ) )

# Make a matrix of Co-efficient of Variation (CV) values
# NOTE: Unlike the stdev calculated above, this is not completely acccurate below, but
# a good enough approximation
CV.matr <- SD.matr / St.matr
# Average fecundity CV of all non-zero F values
res.summ$fec.cv.avg <- mean( CV.matr[ which( fec.matr > 0 ) ] )
# Average survival CV of all non-zero S values
res.summ$surv.cv.avg <- mean( CV.matr[ which( surv.matr > 0 ) ] )
# Average Coefficient of Variation (CV), averaged over all non-zero matrix elements
res.summ$cv.avg <- mean( CV.matr[ which( St.matr > 0 ) ] )

# ----------------------------------------------------------------------------------------------- #
# Variation based on comparing abundance meterics at time step 10
# First check that there are at least 10 time steps
if ( mp.in$MaxDur > 10 ){
  Nmax.10 <- mp.res$PopAll$Max[10]
  Nmin.10 <- mp.res$PopAll$Min[10]
  # Nmax/Nmin for average abundance for time step 10
  res.summ$N.maxVmin.10 <- Nmax.10 / Nmin.10
  # N(+1S.D.)/N(-1S.D.) for average abundance for time step 10
  Nmean.10 <- mp.res$PopAll$Mean[10]
  N10.plus.SD <- Nmean.10 + mp.res$PopAll$StDev[10]
  N10.minus.SD <- Nmean.10 - mp.res$PopAll$StDev[10]
  res.summ$N.plusVmin.SD <- N10.plus.SD / N10.minus.SD
  # CV of N (calculated as (N(+1S.D) - N(-1S.D.)/(2*Navg), from Average Abundance 
  # for time step 10)
  res.summ$N.CV.10 <- (N10.plus.SD - N10.minus.SD)/(2*Nmean.10)  
} else {
  # If there are not at least 10 time steps, then set these values to 'NA'
  res.summ$N.maxVmin.10 <- 'NA'
  res.summ$N.plusVmin.SD <- 'NA'
  res.summ$N.CV.10 <- 'NA'
}
# ----------------------------------------------------------------------------------------------- #
# Determine age/stage of first reproduction
#
# The method used here is a bit complicated, but designed to accomidate
# matrices in which not all fecundities are necessarily in the first row
# of the stage matrix
# Calculate column sums of fecundity matrix
fec.vec <- colSums(fec.matr)
# The first non-zero element of this vector corresponds to the first column with 
# a fecundity value in the stage matrix, which corresponds to the first age/stage
# of reproduction
first.stage <- which( fec.vec > 0 )
first.stage <- first.stage[1]
res.summ$St.First.Rep <- first.stage
# Name of age/stage of first reproduction
res.summ$St.First.Rep.Name <- mp.in$StProp[[first.stage]]$StName

##### BEGIN ENDPOINT CALCULATIONS (Dependent Variables) #####
# ----------------------------------------------------------------------------------------------- #
# Extinction Risk: Risk of total extinction by the final time step 
#
# Get the number of replications that resulted in 'extinction'
ext.rep <- length( which( mp.res$Replications$Ter == 0 ) )
# Exinction risk = (# reps resulting in extinction) / (total # of reps)
res.summ$ext.risk <- ext.rep / mp.res$SimRep # SimRep is the number of replicaitons in this simulation
# ----------------------------------------------------------------------------------------------- #
# Threshold: The quasi-extinction threshold. Value set in Stochasticity section of MP Module.
# Used to calculate time to quasi-extinction.
res.summ$threshold <- mp.in$ExtinctThr
# ----------------------------------------------------------------------------------------------- #
# Prob(maxt): Risk of falling below the Threshold at least once by the final time step
# Calculated using information regarding the time of first crossing of the extinction threshold
res.summ$prob.thresh.maxt <- sum( mp.res$TimeCross$QuasiExtinct ) / mp.res$SimRep
# ----------------------------------------------------------------------------------------------- #
# MedTime: Median time to fall below the quasi-extinction threshold
#
# First calculate cumulative sum of TimeCross$QuasiExtinct / SimRep (i.e. the prob of first 
# crossing the threshold at that time step)
prob.cross <- cumsum( mp.res$TimeCross$QuasiExtinct / mp.res$SimRep )
if ( any( prob.cross >= 0.5 ) ) { # Check if any values are greater than 0.5
  # Get the minimum time-step for which the prob.cross is greater than or equal to 0.5. In most cases
  # this will be an overestimate of the median time to quasi-extinction
  #browser()
  med.time.temp <- min( which(prob.cross >= 0.5) )
  # Approximation for med.time if med.time.temp does not occur at exactly 0.5
  if ( prob.cross[med.time.temp] != 0.5 ) {
    if ( med.time.temp == 1) {
      med.time <- 0.5/prob.cross[med.time.temp]
    } else {
      med.time <- (0.5-prob.cross[med.time.temp-1])/(prob.cross[med.time.temp]-prob.cross[med.time.temp-1]) + (med.time.temp-1)
    }
  } else {
    med.time <- med.time.temp
  }
  # Correct for the 'stepsize'
  res.summ$med.time.cross <- med.time * mp.in$stepsize
} else {
  # If the median time to cross does not occur before the end of the simulation duration, 
  # then set to the maximum simulation duration + 1 (as done in RAMAS ResultsOut)
  res.summ$med.time.cross <- mp.in$MaxDur + 1
}
# ----------------------------------------------------------------------------------------------- #
# Prob(50): Risk of falling below a population size of 50 individuals at least once by the final
# time step
res.summ$prob.50 <- length( which( mp.res$Replications$Min <= 50 ) ) / mp.res$SimRep
# ----------------------------------------------------------------------------------------------- #
# Prob(250): Risk of falling below a population size of 250 individuals at least once by the final
# time step
res.summ$prob.250 <- length( which( mp.res$Replications$Min <= 250 ) ) / mp.res$SimRep
# ----------------------------------------------------------------------------------------------- #
# Prob(1000): Risk of falling below a population size of 1000 individuals at least once by the final
# time step
res.summ$prob.1000 <- length( which( mp.res$Replications$Min <= 1000 ) ) / mp.res$SimRep
# ----------------------------------------------------------------------------------------------- #
# Prob(50%): Risk of falling below a population size of 50% the intial metapop size at least once 
# by the final time step
res.summ$prob.Thr.50 <- length( which( mp.res$Replications$Min <= Thr.50 ) ) / mp.res$SimRep
# ----------------------------------------------------------------------------------------------- #
# ExpMinN: Expected minimum total abundance (also known as EMA)
res.summ$exp.min.n <- mean( mp.res$Replications$Min )
# ----------------------------------------------------------------------------------------------- #
# SdErr.ExpMinN: Standard error of EMA
res.summ$sderr.ema <- sqrt( var( mp.res$Replications$Min ) ) / sqrt( mp.res$SimRep )
# ----------------------------------------------------------------------------------------------- #
# N(maxt): Mean metapopulation abundance at the final time step maxt
# This value is the same as mp.res$PopAll$Mean( maxt )
res.summ$n.mean <- mean( mp.res$Replications$Ter )
# ----------------------------------------------------------------------------------------------- #
# SD of N: Standard deviation of metapopulation occupancy at the final time step maxt
# This value is the same as mp.res$PopAll$StDev( maxt )
res.summ$n.stdev <- sqrt( var(mp.res$Replications$Ter) )
# ----------------------------------------------------------------------------------------------- #
# % Metapop Ab Change: Percent change in total metapopulation size compared to initial 
# metapopulation size
res.summ$metapop.chng <- ( res.summ$n.mean - res.summ$metapop.initab ) / res.summ$metapop.initab
# ----------------------------------------------------------------------------------------------- #
# Occ(maxt): Mean metapopulation occupance at the final time step maxt
res.summ$occ.maxt <- mp.res$Occupancy$Mean[ mp.in$MaxDur ]
# ----------------------------------------------------------------------------------------------- #
# SD of Occ: Standard deviation of metapopulation occupancy at the final time step maxt
res.summ$occ.stdev.maxt <- mp.res$Occupancy$StDev[ mp.in$MaxDur ]
# ----------------------------------------------------------------------------------------------- #
# Quantiles of Terminal Population Trajectory.  5%, 25%, 50%, 75%, and 95% quantiles
quants <- quantile( mp.res$Replications$Ter, prob = c(0.05, 0.25, 0.50, 0.75, 0.95) )
res.summ$quant.05 <- quants[1]
res.summ$quant.25 <- quants[2]
res.summ$quant.50 <- quants[3]
res.summ$quant.75 <- quants[4]
res.summ$quant.95 <- quants[5]
# ----------------------------------------------------------------------------------------------- #
# Harvest Results: Average total harvest, St.Dev. of total harvest, Minimum of total harvest
# Maximum of total harvest
res.summ$harv.avg <- mp.res$HarvestTot$Mean
res.summ$harv.stdev <- mp.res$HarvestTot$StDev
res.summ$harv.min <- mp.res$HarvestTot$Min
res.summ$harv.max <- mp.res$HarvestTot$Max

# ----------------------------------------------------------------------------------------------- #
###################################################################################################
#### Begin section concering Spatial Endpoints ####
# NOTE: Incorrect indentation for this 'if' statement
if ( spatial ) {
print('mp.results: Begin reading SPATIAL endpoints')
# ----------------------------------------------------------------------------------------------- #
# Dispersal Measures from Dispersal Distance Function
# 
# Use Dispersal Distance Function?
res.summ$use.disp.dist.func <- mp.in$UseDispDistFunc
# Get 'a' value from Disp. Dist. Func.
res.summ$disp.dist.func.a <- mp.in$DispDistFunc[1]
# Get average dispersal distance, b value from Disp. Dist. Func.
res.summ$avg.disp.dist.b <- mp.in$DispDistFunc[2]
# Get 'c' value from Disp. Dist. Func
res.summ$disp.dist.func.c <- mp.in$DispDistFunc[3]
# Get max dispersal distance, Dmax value from Disp. Dist. Func.
res.summ$max.disp.dist.Dmax <- mp.in$DispDistFunc[4]
# ----------------------------------------------------------------------------------------------- #
# Dispersal Measures from Dispersal Matrix
#
# Define a new variable that is the dispersal matrix to make code clearer
disp.matr <- mp.in$DispMatr
# Average total dispersal rate (dispersal matrix: sum of each column, averaged over columns)
#
# We will consider dispersal for only those populations that exist at time-step zero 
# Determine pops that exist at time-step zero
pops.t0 <- which( mp.in$PopData_df$K > 0 )
# Define this smaller dispersal matrix
disp.matr.sub <- disp.matr[pops.t0,pops.t0]
# Take the sum of each column, this is essentially the percent of individuals that leave a 
# given population
disp.t0.colSums <- apply( disp.matr.sub, 2, sum)
# Mean dispersal (from a population) rate
res.summ$mean.t0.disp.rate <- mean( disp.t0.colSums )
# Lower and Upper quartile of total dispersal rate
disp.t0.rate.quant <- quantile( disp.t0.colSums, prob = c(0.25,0.75) )
res.summ$loquart.t0.disp.rate <- disp.t0.rate.quant[1]
res.summ$hiquart.t0.disp.rate <- disp.t0.rate.quant[2]

# We will now consider dispersal for only those populations with non-zero initial abundance
# and  that exist at time-step zero 
# Determine pops with initial abundance greater than 0 at time-step zero
pops.t0.Ngt0 <- which( mp.in$PopData_df$InitAbund > 0 )
# Mean dispersal (from a population) rate
res.summ$mean.t0.Ngt0.disp.rate <- mean( disp.t0.colSums[ pops.t0.Ngt0 ] )
# Lower and Upper quartile of total dispersal rate
disp.t0.Ngt0.rate.quant <- quantile( disp.t0.colSums[ pops.t0.Ngt0 ], prob = c(0.25,0.75) )
res.summ$loquart.t0.Ngt0.disp.rate <- disp.t0.Ngt0.rate.quant[1]
res.summ$hiquart.t0.Ngt0.disp.rate <- disp.t0.Ngt0.rate.quant[2]

# Average dispersal rate to nearest neighbor (assume that max dispersal in each column is
# to nearest neighbor)
disp.matr.colMax <- apply( disp.matr.sub, 2, max )
res.summ$mean.t0.nearest.disp.rate <- mean( disp.matr.colMax )
# Lower and Upper quartile of nearest neighbor dispersal rate
nearest.disp.rate.quant <- quantile( disp.matr.colMax, prob = c(0.25,0.75) )
res.summ$loquart.t0.nearest.disp.rate <- nearest.disp.rate.quant[1]
res.summ$hiquart.t0.nearest.disp.rate <- nearest.disp.rate.quant[2]

# Average dispersal rate to nearest neighbor, conditioned on a non-zero
# initial population abundance
res.summ$mean.t0.Ngt0.nearest.disp.rate <- mean( disp.matr.colMax[ pops.t0.Ngt0 ] )
# Lower and Upper quartile of nearest neighbor dispersal rate
nearest.disp.rate.Ngt0.quant <- quantile( disp.matr.colMax[ pops.t0.Ngt0 ], prob = c(0.25,0.75) )
res.summ$loquart.t0.Ngt0.nearest.disp.rate <- nearest.disp.rate.Ngt0.quant[1]
res.summ$hiquart.t0.Ngt0.nearest.disp.rate <- nearest.disp.rate.Ngt0.quant[2]


# ----------------------------------------------------------------------------------------------- #
# Correlation Measures
#
# Use Correlation Distance Function?
res.summ$use.corr.dist.func <- mp.in$UseCorrDistFunc
# Average correlation distance, based on b value of corr. dist. func.
res.summ$avg.corr.dist.b <- mp.in$CorrDistFunc[2]

# ----------------------------------------------------------------------------------------------- #

} # End 'if (spatial)' statement

# ----------------------------------------------------------------------------------------------- #
###################################################################################################
#### Begin section of endpoints from PTC file(s) ####
## Here the user can specify ptcFiles as being one file or a concatenated vector of several
## *.ptc files.
## All of the below indices will be extracted for patches with non-zero abundance
# NOTE: Incorrect indentation for this 'if' statement
if ( spatial & ptc ) {

# Inform user that endpoints from PTC file are being extracted
print( paste('mp.results: Begin extraction of endpoints from *.ptc file: ', ptcFiles))

## Check that the ptcFiles exist
ptcFilesExist <- lapply(ptcFiles, file.exists)
ptcFilesExist <- unlist(ptcFilesExist)
# If any file does not exist, then stop the program
if( !all(ptcFilesExist) ){
  stop('Did not find all *.ptc file(s).')
}

## Check that there are an equal number of ptcFiles and ptcFileIter
if( !(length(ptcFiles)==length(ptcFileIter)) ) {
  stop('Uneven number of ptcFiles and ptcFileIter.')
}

## Read the ptcFile(s)
ptc.files <- lapply( ptcFiles, ptc.read )

## Define the column names for a table of Isolation and Area and Shape Indices. These
## column names are used in both habdyn = TRUE or FALSE conditions
isolation.colnames <- c('metapop.abundance', 'avg.dist.occupied.pops',
  'avg.dist.near.neighbor', 'low.quart.dist.near.neighbor',
  'upper.quart.dist.near.neighbor', 'max.dist.near.neighbor',
  'frag.index.hra')
area.shape.colnames <- c('patch.n', 'tot.patch.area', 'tot.core.area',
  'tot.edge.length','lg.patch.area','lg.patch.core','overall.edge.area.ratio',
  'overall.shape.ind','overall.frac.dim','avg.patch.area','avg.core.area',
  'avg.edge.length','avg.shape.ind','avg.frac.dim')
indice.colnames <- c(isolation.colnames, area.shape.colnames)

## Create two empty vectors in which data will be added to based on the number 
## *.ptc files (i.e. iterations of the ptc.cnt for loop)
# Indice column names vector
ind.coln.full <- vector()
# Indices vector
ind.full <- vector()

## Begin a for loop that goes through each of the *.ptc files included
for ( ptc.cnt in 1:length(ptcFiles) ){
  # Look at the ptc file and determine the number of patches that exist
  # in the current file
  PatchN <- ptc.files[[ ptc.cnt ]]$PopN
  
  if (habdyn) {
    ## First, let's consider the scenario in which habdyn = TRUE
    ## In this scenario the HabDyn Module in RAMAS GIS was used to incorporate 
    ## habitat dynamics into the *.mp file.

    ## Check that hdhFile exists
    if( !file.exists( hdhFile ) ){
      stop( paste( 'Did not find HabDyn History file: ', hdhFile ) )
    }  
    ## Read the HabDyn History file
    hdhist <-hdhist.read( hdhFile )
    
    ## Special section for NASA Project! EDIT or Remove for Generic!!!
    lastIter <- max(hdhist$patch.info.mat$iter)

    ## Save the original 'ptcFileIter' vector
    ptcFileIter.orig <- ptcFileIter
    
    if( ptcFileIter[ptc.cnt] > lastIter) { ptcFileIter[ptc.cnt]<-lastIter }
    
    # Look at hdhist file and determine which *.mp populations exist at 
    # ptcFileIter[ptc.cnt] iteration (time step).  There can be a greater number of
    # mp populations than ptc populations, since some populations in the mp file will
    # have K=0 for some time steps
    patch.mat.sub <- subset( hdhist$patch.info.mat, 
      subset=hdhist$patch.info.mat$iter==ptcFileIter[ptc.cnt])
    # Limit this subset to only the patches that exist in the ptc file for this 
    # iteration step
    patch.mat.sub <- patch.mat.sub[1:PatchN,]
    # Address the case where there is a zero in the new2old column (usually 
    # only for the first iteration step)
    iter.mp.pops <- ifelse( test=patch.mat.sub$new2old==0,
      yes=patch.mat.sub$patch, no=patch.mat.sub$new2old)
  
    # Determine the matching mp simulation start year for a particular iteration
    iter.start.yr <- unique( patch.mat.sub$Fyear ) 
    if ( length(iter.start.yr) > 1 ){
      stop('mp.results: HabDyn History file is correctly formated.
        There should only be one Fyear per iter step.')
    } 
    
    ### WARNING: NASA Specific addition
    ## Assume that the iter.start.yr corresponds to the ptcFileIter.orig value
    ## minus 1.  The minus 1 comes from the fact that time-steps start at 1 and
    ## years start from 0.
    iter.start.yr <- ptcFileIter.orig[ ptc.cnt ] - 1
    
    # Make a new set of col names that includes the start year  for this ptc file
    indice.colnames.temp <- cbind( indice.colnames, rep( ptcFileIter[ptc.cnt],
      length( indice.colnames)))
    indice.colnames.temp <- apply( indice.colnames.temp,1,paste,collapse='.')
        
    # Now check which populations have non-zero abundance.  If iter.start.yr==0, then use
    # initial abundance from PopData_df.  If iter.start.yr>1, then use Mean Pop Trajectory value
    # for the individual populations from the mp$results$PopInd section
    
    ###browser()
    
    if ( iter.start.yr==0 ){
      pop.ab <- mp.in$PopData_df[ iter.mp.pops, 4 ] # InitAb is column 4
    } else {
      pop.ab <- mp.res$PopInd[ iter.start.yr, 1, iter.mp.pops ]
    }
    # Calculate the total metapop.abundance for this 'iter' time step.
    metapop.abundance <- sum(pop.ab)
        
    # Determine the ptc populations with non-zero abundance
    patch.nonzero <- patch.mat.sub$patch[ which( pop.ab > 0 ) ]
    
  } else {
    ## If habdyn=FALSE, then it is assumed the that there is only one *.ptc file associated
    ## with the *.mp file. That is, no habitat dynamics were calculated using the HabDyn
    ## module. Using HabDyn is NOT the only way to incorporate habitat dynamics (e.g. simple
    ## declines in K can be added to the *.mp file directly), but for the sake of this 
    ## calculation, we are assuming that if there is no HabDyn History file included, then
    ## the *.ptc file is associate with only year 0 (initial year) of the simulation.
    
    # Assume *.ptc file corresponds with year 0 of the *.mp model
    iter.start.yr = 0
    
    # Make a new set of col names that includes the start year  for this ptc file
    indice.colnames.temp <- cbind( indice.colnames, rep( ptcFileIter[ptc.cnt],
      length( indice.colnames)))
    indice.colnames.temp <- apply( indice.colnames.temp,1,paste,collapse='.')
    
####    # Determine the ptc populations with non-zero abundance
####    # In the scenario where HabDyn = FALSE, the population numbers in the *.ptc
####    # and *.mp files should match.  Additionally, in this case we are relying on 
####    # initial population abundance.
####    pop.ab <- mp.in$PopData_df[ 1:PatchN, 4 ] # InitAb is column 4
####    patch.nonzero <- which( pop.ab > 0 )

#### WARNING: Special changes here just for the NASA Project ####
#### These changes assume that there is a vector of ptc files and
#### ptcFileIter numbers.  These changes allow us to use the same ptc
#### file, namely the first time-step ptc file, to calculate the spatial
#### structure for several different time-steps, in order to facilitate the
#### comparisons used in the analysis of the NASA Project.
#### Uncomment lines above for generic version of this script
    if ( ptcFileIter[ptc.cnt] == 1 ) {
      # If iter is 1, use initial abundance
      pop.ab <- mp.in$PopData_df[ 1:PatchN, 4 ] # InitAb is column 4
    } else {
      # Use the time step pop abundance from the mp file
      # in this case, the iteration is the year + 1, so iter 21
      # corresponds to mp results from year 20
      pop.ab <- mp.res$PopInd[ (ptcFileIter[ptc.cnt]-1), 1, 1:PatchN ]
    }
    metapop.abundance <- sum(pop.ab)

    patch.nonzero <- which( pop.ab > 0 )
  } #END IF HABDYN
  
  ## Account for situations in which patch.nonzero = 0
  if ( length(patch.nonzero) == 0 ){
    isolation.indices.iter <- c(0,rep(NA,6))

    ## Area and Shape Indices vector for this HabDyn iteration
    area.shape.indices.iter <- c( 0, 0, 
                                  0, 0, 0, 0,
                                  NA, NA, NA, 
                                  0, 0, 0, NA,
                                  NA)
    indices.iter <- c(isolation.indices.iter, area.shape.indices.iter)
    
  } else {    
    ## Calculate Isolation Indices based on distance matrix of non-zero abundance
    ## populations
    #
    if ( length(patch.nonzero) > 1){
      # Calculate distance matrix for non-zero populations
      dist.matr.nonzero <- ptc.files[[ ptc.cnt ]]$DistMatr[patch.nonzero,patch.nonzero]
      # Convert 0 to NA
      dist.matr.nonzero[ dist.matr.nonzero==0 ] <- NA
      
      # Calculate minimum values for each column
      dist.min.nonzero <- apply( dist.matr.nonzero, 2, min, na.rm=TRUE )      
      # Calculate mean values for each column
      dist.mean.nonzero <- apply( dist.matr.nonzero, 2, mean, na.rm=TRUE )

      # Average distance to other occupied populations
      avg.dist.occupied.pops <- mean( dist.mean.nonzero )     
      # Average distance to nearest neighbor
      avg.dist.near.neighbor <- mean(dist.min.nonzero)
      
      # Lower and Upper quartiles of distance to nearest neighbor
      quartiles.dist.near.neighbor <- quantile( dist.min.nonzero, prob=c(0.25,0.75) )
      lo.quart.dist.near.neighbor <- quartiles.dist.near.neighbor[1]
      hi.quart.dist.near.neighbor <- quartiles.dist.near.neighbor[2]
      
      # Maximum distance to nearest neighbor 
      max.dist.near.neighbor <- max( dist.min.nonzero )
      
      # ----------------------------------------------------------------------------------------------- #
      ## Fragmentation index developed by HRA - 
      ## Percent of total metapopulation occupying populations with immigration rate <= 5%
      ## and abundance <= 50 individuals
      
      # First convert distance matrix to dispersal matrix
      disp.matr.nonzero <- mp.in$DispDistFunc[1] * exp( -(dist.matr.nonzero^mp.in$DispDistFunc[3])/mp.in$DispDistFunc[2] )
      disp.matr.nonzero[ dist.matr.nonzero > mp.in$DispDistFunc[4] ] <- 0
      
      # Extract Population Abundance values for non-zero populations
      pop.ab.nonzero <- pop.ab[ patch.nonzero ]
      # Create matrix of the above vector
      pop.ab.nonzero.mat <- matrix( rep( pop.ab.nonzero, length(patch.nonzero) ), nrow=length(patch.nonzero), byrow=TRUE )
      
      # Calculate the number of dispersers from/to each population
      dispersers <- disp.matr.nonzero * pop.ab.nonzero.mat
      
      # Calculate the number of IMMIGRANTS for each population (row sum of the dispersers matrix)
      immigrants <- apply( dispersers,MARGIN=1,FUN=sum,na.rm=TRUE)
      
      # Calculate the immigration rate
      immigrant.rate <- immigrants/pop.ab.nonzero
      
      # Calculate fragmented patches
      frag.patch <- which( immigrant.rate < 0.05 & pop.ab.nonzero < 50 )
      
      # Calculate population abundance within fragmented patches
      pop.ab.frag.patch <- sum( pop.ab.nonzero[ frag.patch ] )

      # Calculate percent of population in fragmented patches versus total population
      frag.index.hra <- pop.ab.frag.patch / metapop.abundance
            
      #browser() ### DEBUG LINE ###
      # ----------------------------------------------------------------------------------------------- #
      
      ## Isolation Indices Vector for this HabDyn iteration
      isolation.indices.iter <- c( metapop.abundance, avg.dist.occupied.pops,
                                   avg.dist.near.neighbor, 
                                   lo.quart.dist.near.neighbor,
                                   hi.quart.dist.near.neighbor, max.dist.near.neighbor,
                                   frag.index.hra)
    } else {
      ## If there is only ONE non-zero patch, set isolation indices to default values

      ## If patch.nonzero = 1, then frag.index.hra should be returned as either 1 or 
      ## 0, depending on the population size.  If the population size is small, then 
      ## return 1 (completely fragmented), if large then return 0.
      frag.index.hra <- ifelse( metapop.abundance < 50, 1, 0 )
      
      isolation.indices.iter <- c(metapop.abundance, rep(NA,5), frag.index.hra)
    }
    
    
    ## Calculate Area and Shape Indices for non-zero abundance populations
    #
    # Number of non-zero patches
    patch.n.nonzero <- length( patch.nonzero )
    
    # Extract non-zero population (patch) information from the 
    # ptc.file$PopLandInd_df data frame and from ptc.file$PopPatchChar_df
    PopLandInd_nonzero <- ptc.files[[ptc.cnt]]$PopLandInd_df[ patch.nonzero, ]
    PopPatchChar_nonzero <- ptc.files[[ptc.cnt]]$PopPatchChar_df[ patch.nonzero, ]
    
    # Total area of patches (km^2)
    # NOTE: must convert from area in cells to km^2 by multiplying the number
    # of cells by the cell length squared
    tot.patch.area.cells <- sum(PopPatchChar_nonzero$Area)
    tot.patch.area.km2 <- tot.patch.area.cells * ptc.files[[ptc.cnt]]$cell.length^2
    
    # Total core area (km^2)
    tot.core.area.cells <- sum(PopLandInd_nonzero$core.area)
    tot.core.area.km2 <- tot.core.area.cells * ptc.files[[ptc.cnt]]$cell.length^2
    
    # Total edge length (km)
    tot.edge.length <- sum(PopLandInd_nonzero$perimeter)
    
    # Area of largest patch
    lg.patch.area.cells <- max(PopPatchChar_nonzero$Area)
    lg.patch.area.km2 <- lg.patch.area.cells * ptc.files[[ptc.cnt]]$cell.length^2
    
    # Core area of the largest patch
    # Identifies largest core of largest patch by finding largest patch, just in case 
    # largest core is not the same as largest patch
    lg.patch.core.cells <- PopLandInd_nonzero$core.area[ which( PopPatchChar_nonzero$Area == lg.patch.area.cells ) ]
    lg.patch.core.km2 <- lg.patch.core.cells * ptc.files[[ptc.cnt]]$cell.length^2
    
    # Overall edge:area ratio (1/km)
    # Total Edge Length / Total Area of Patches
    overall.edge.area.ratio <- tot.edge.length / tot.patch.area.km2
    
    # Overall shape index
    # (0.25 * Total Edge Length) / Sqrt( Total Area )
    overall.shape.ind <- (0.25*tot.edge.length)/sqrt(tot.patch.area.km2)
    
    # Overall fractal dimension
    # 2 * ln( 0.25*Total Edge Length) / ln( Total Area)
    overall.frac.dim <- 2 * log( 0.25*tot.edge.length ) / log(tot.patch.area.km2)
    
    # Average patch area
    avg.patch.area.cells <- mean( PopPatchChar_nonzero$Area )
    avg.patch.area.km2 <- avg.patch.area.cells * ptc.files[[ptc.cnt]]$cell.length^2
    
    # Average core area
    avg.core.area.cells <- mean( PopLandInd_nonzero$core.area )
    avg.core.area.km2 <- avg.core.area.cells * ptc.files[[ptc.cnt]]$cell.length^2
    
    # Average edge length
    avg.edge.length <- mean( PopLandInd_nonzero$perimeter )
    
    # Average shape index
    avg.shape.ind <- mean( PopLandInd_nonzero$shp.index )
    
    # Average fractal dimension
    avg.frac.dim <- mean( PopLandInd_nonzero$fractal.dim )
    
    ## Area and Shape Indices vector for this HabDyn iteration
    area.shape.indices.iter <- c( patch.n.nonzero, tot.patch.area.km2, 
                                  tot.core.area.km2, tot.edge.length, lg.patch.area.km2, lg.patch.core.km2,
                                  overall.edge.area.ratio, overall.shape.ind, overall.frac.dim, 
                                  avg.patch.area.km2, avg.core.area.km2, avg.edge.length, avg.shape.ind,
                                  avg.frac.dim)
    
    ## Concatenate isolation indices and area/shape indices
    indices.iter <- c( isolation.indices.iter, area.shape.indices.iter )
    
  } # End if patch.nonzero == 0 conditional 
  
  
  ## Concatenate full indice vectors
  ind.coln.full <- c( ind.coln.full, indice.colnames.temp )
  ind.full <- c( ind.full, indices.iter)
  
  
  ptc.file.temp <- ptc.files[[ ptc.cnt ]]

} # End for ptc.cnt loop

# Convert full lists into a data frame
patch.indices.df <- as.data.frame( t(ind.full) )
colnames( patch.indices.df ) <- ind.coln.full
# Convert this data frame into a list and add to res.summ
patch.indices.list <- as.list( patch.indices.df )
res.summ <- c( res.summ, patch.indices.list )
###browser()

} # End 'if (spatial)' statement

###res.summ.df <- as.data.frame(t(res.summ))
res.summ.df <- as.data.frame(res.summ)
row.names(res.summ.df) <- mpFile
return(res.summ.df)
}
