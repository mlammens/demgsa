ptc.read <- function( ptcFile ) {
  ########################################################################################
  ## Read data from a *.ptc file, used in RAMAS GIS - 
  ## Spatial Data module
  ##
  ## Args:
  ##  ptcFile: The path and name of a *.ptc file to be read.
  ##
  ########################################################################################
  
  # Print the name of the ptc file to the console
  print( paste( "Begin ptc.read function with file: ", ptcFile ) )
  # Read the *.ptc file in as a long unsorted 'list' structure
  Ptc <- readLines( ptcFile )
  # Clear first line of "map=\xff". This step is necessary to avoid errors when using 
  # regex functions in this script
  Ptc[1] <- sub("(map=).*","\\1",Ptc[1])

  # Set values that will be used through out the function
  MigrationLine <- which(Ptc == "Migration")
  CorrLine <- which(Ptc == "Correlation")
  
  # Check if Results are present.  If yes, set ResLine and LandLine variables
  if ( any(grepl( "^Results", Ptc )) ) {
    print( paste("ptc.read: Results present in *.ptc file: ", ptcFile) )
    ResLine <- grep("Results",Ptc)
    LandLine <- grep("Landscape indices",Ptc)
  } else {
    print( paste("ptc.read: No Results found in *.ptc file: ", ptcFile) )
  }

  # Create an empty 'list' structure to store ptc file information
  Ptc.list <- vector("list",length=0)
  
  ## As of RAMAS GIS 5.1, the first 28 lines of the *.ptc file are fixed.  The current 
  ## version of this program assumes that within these first 28 lines to location (i.e.
  ## line number) of some information is fixed.
  #
  # Landscape input file version: Check that file is a 'Landscape' file and, if yes, get 
  # first two numbers in string that are separated by a period using regexs
  if ( grepl("^Landscape", Ptc[1]) ) {
    Ptc.list$version <- sub(".*([0-9]).([0-9]).*","\\1\\2",Ptc[1])
  } else {
    stop("Error reading *.pct file. Check that this is a correct *.ptc file")
  }
  # PTC Scenario Title
  Ptc.title <- Ptc[2]
  # Comment Lines: Four lines of comments in *.ptc file
  Ptc.list$comments <- Ptc[3:6]
  # Cell length
  Ptc.list$cell.length <- as.numeric( Ptc[7] )
  # Habitat Suitability Function
  Ptc.list$HSI.func <- Ptc[8]
  
  ## Habitat Relationships (HSI) section
  # Habitat suitability threshold for patches
  Ptc.list$HSI.Threshold <- as.numeric( Ptc[10] )
  # Neighborhood Distance (cells)
  Ptc.list$Neigh.Dist <- as.numeric( Ptc[11] )
  # Habitat suitability map color and Input is patch map boolean
  color.pmap <- unlist( strsplit( Ptc[12], "," ) )
  Ptc.list$HSI.MapColor <- color.pmap[1]
  Ptc.list$HSI.IsPatchMap <- as.logical( color.pmap[2] )
  # Number of decimals to export in habitat suitability map
  Ptc.list$HSI.ExportDecN <- as.numeric( Ptc[13] )
  
  ## Link to Metapopulation (MP) section
  # Carrying Capacity (K) function
  Ptc.list$MP.K.func <- Ptc[14]
  # Rmax function
  Ptc.list$MP.Rmax.func <- Ptc[15]
  # Initial abundance function
  Ptc.list$MP.InitN.func <- Ptc[16]
  # Relative fecundity function
  Ptc.list$MP.RelFec.func <- Ptc[17]
  # Relative survival function
  Ptc.list$MP.RelSur.func <- Ptc[18]
  # *.mp file to get other data from
  Ptc.list$MP.DefaultMP <- Ptc[20]
  # Catastrophe 1 - Local Prob
  Ptc.list$MP.Cat1.LocProb <- Ptc[21]
  # Catastrophe 1 - Local Mult
  Ptc.list$MP.Cat1.LocMult <- Ptc[22]
  # Catastrophe 2 - Local Prob
  Ptc.list$MP.Cat2.LocProb <- Ptc[23]
  # Catastrophe 2 - Local Mult
  Ptc.list$MP.Cat2.LocMult <- Ptc[24]
  # Habitat-based distances (from a friction map) (Yes/No)
  Ptc.list$MP.HabitatBasedDistances <- Ptc[25]
  # Friction map file name
  Ptc.list$MP.FrictionMapFilename <- Ptc[26]
  # Distance calculation (Edge to Edge, Center to Edge, or Center to Center)
  Ptc.list$MP.DistCalc <- Ptc[27]
  
  ## Input Maps Information
  # Number of input maps
  Ptc.list$MapN <- as.numeric(Ptc[28])
  ## For each map, read in information
  # Create an empty AllMapData list
  AllMapData <- vector('list',length=0)
  if ( Ptc.list$MapN > 0 ){
    # Set the number of lines of map information.  If the number of maps is non-zero
    # then map info will begin at line 29 and continue to 28+(MapN*5)
    MapLines <- Ptc[29:(28+(Ptc.list$MapN*5))]
    # (Re)Create an empty MapData list
    MapData <- vector('list',length=0)
    for ( map in 1:Ptc.list$MapN ) {
      LineReadOffset <- (Ptc.list$MapN-1)*5
      MapData$Name <- MapLines[ 1 + LineReadOffset ]
      MapData$FileName <- MapLines[ 2 + LineReadOffset ]
      MapData$FileFormat <- MapLines[ 3 + LineReadOffset ]
      MapData$Color <- MapLines[ 4 + LineReadOffset ]
      MapData$ColN <- as.numeric( MapLines[ 5 + LineReadOffset ] )
      # Add new MapData list to AllMapDataList
      AllMapData[[ map ]] <- MapData
    }    
  } else {
    print( 'No maps included in this *.ptc file.')
  }
  Ptc.list$MapData <- AllMapData
  
  ## Default Population Information
  # Note default population information, assuming it is the line immediately preceding
  # the MigrationLine.
  Ptc.list$DefaultPop <- Ptc[ MigrationLine - 1 ]
  
  ## Migration (Dispersal) and Correlation Distance Functions
  # UseDispDistFunc: True if dispersal rates are based on dispersal distance function; false if
  # they are specified in the dispersal matrix. NOTE: For *.ptc files, this value should
  # always be TRUE
  Ptc.list$UseDispDistFunc <- as.logical( Ptc[MigrationLine + 1] )
  # DispDistFunc: Dispersal-distance function parameters - a, b, c, Dmax - Mij = a exp(-Dij^c/b)
  Ptc.list$DispDistFunc <- as.numeric(unlist(strsplit( Ptc[MigrationLine + 2],',' )))
  # UseCorrDistFunc: True if correlations between populations is based on correlation distance
  # function. NOTE: For *.ptc files, this value should always be TRUE
  Ptc.list$UseCorrDistFunc <- as.logical( Ptc[CorrLine +1] )
  # CorrDistFunc: Correlation-distance function parameters - a, b, c - Cij = a exp(-Dij^c/b)
  Ptc.list$CorrDistFunc <- as.numeric(unlist(strsplit( Ptc[CorrLine + 2],',' )))
  
  ##------------------------------------------------------------------------------------##
  ## Below are data associated with the Results section of a *.ptc file, and while only
  ## be present in those *.ptc files that have been run.  This section includes information
  ## on the HSI Histogram, which is actually noted prior to the Results line.
  
  if ( any(grepl( "^Results", Ptc )) ) {  
    
    ## HSI Histogram Information:
    # Determine the HSI map line
    hsi.map.line <- grep("^HSI map",Ptc)
    ### FUTURE DEV - Make this information useful in R
    Ptc.list$HSI.Histogram <- Ptc[ (hsi.map.line+1):(ResLine-2) ]
    Ptc.list$HSI.MapInfo <- Ptc[ (ResLine-1) ]
    
    # Results Date Information
    Ptc.list$ResultsDate <- Ptc[ ResLine ]
    
    # Number of populations
    PopLine <- unlist( strsplit( Ptc[ ResLine + 1 ]," " ) )
    Ptc.list$PopN <- as.numeric(PopLine[1])
    
    ## Population Information: This section is similar to the format employed in 
    ## RAMAS Version 5.0.  
    # Define Population Data Row Names
    PopData_df_colnames <- c("name","X_coord","Y_coord","InitAbund","DensDep","MaxR","K","Ksdstr","Allee","KchangeSt","DD_Migr","Cat1.Multiplier","Cat1.Prob","IncludeInSum","StageMatType","RelFec","RelSur","localthr","Cat2.Multiplier","Cat2.Prob","SDMatType","TargetPopK","Cat1.TimeSinceLast","Cat2.TimeSinceLast","RelDisp")
    # Read Population Data as data.frame
    PopData_df <- read.csv( ptcFile, header=FALSE, skip = (ResLine+1), nrows=Ptc.list$PopN )
    PopData_df <- PopData_df[1:length(PopData_df_colnames)]
    # Turn NAs into empty strings
    PopData_df[is.na(PopData_df)] <- ''
    # Assign columns names
    names(PopData_df) <- PopData_df_colnames
    Ptc.list$PopData_df <- PopData_df
    
    ## Population Patch Characteristics
    ## Total Habitat Suitability (HSI), Patch Area (#Cells), pminX, pminY, pmaxX, pmaxY,
    ## Population (patch) ID
    # PopPatchChar Row Names
    PopPatchChar_df_colnames <- c("HSI","Area","pminX","pminY","pmaxX","pmaxY","PopID")
    PopPatchChar_df <- read.table( ptcFile, header=FALSE, skip = (ResLine+Ptc.list$PopN+3), nrows=Ptc.list$PopN )
    names(PopPatchChar_df) <- PopPatchChar_df_colnames
    Ptc.list$PopPatchChar_df <- PopPatchChar_df
    
    ## Summary Line
    ## Not sure what to do with this line yet.  It includes the following values:
    ## Cell length; Max Dispersal Distance; Max Patch Area; Unidentified 1; Unidentified 2
    Ptc.list$SummaryLine <- Ptc[ ResLine + (2*Ptc.list$PopN) + 4 ]
    
    ## Distance Matrix
    ## If center-to-center or edge-to-edge distance is selected, then the matrix in the 
    ## *.ptc file is a lower triangular matrix, and reading it in is more complicated 
    ## than simply reading in a square matrix.
    # If matrix is square (i.e. "Center to edge")
    if ( Ptc.list$MP.DistCalc == "Center to edge" ) {
      DistMatr <- read.table( ptcFile, header=FALSE, skip = (ResLine + (2*Ptc.list$PopN) + 4), nrows=Ptc.list$PopN)
      DistMatr <- as.matrix(DistMatr)
      Ptc.list$DistMatr <- DistMatr
    } else if ( any( Ptc.list$MP.DistCalc == c("Edge to edge", "Default: Center to center"))) {
      # Use fill function, assuming there are PopN columns
      DistMatr <- read.table( ptcFile, header=FALSE, skip = (ResLine + (2*Ptc.list$PopN) + 4), nrows=Ptc.list$PopN, fill=TRUE, col.names = (1:Ptc.list$PopN))
      DistMatr <- as.matrix(DistMatr)
      DistMatr[is.na(DistMatr)] <- 0
      DistMatr <- DistMatr + t(DistMatr)
      Ptc.list$DistMatr <- DistMatr
    }
    
    ## Landscape Indices - Populations
    ## Five indices for each population (patch):
    ## Perimeter, Shape Index, Fractal Dim., Core Area, Edge-to-Area Ratio
    # PopLandInd Row Names
    PopLandInd_df_colnames <- c("perimeter","shp.index","fractal.dim","core.area","edge.area.ratio")
    Ptc.list$PopLandInd_df <- read.table( ptcFile, header=FALSE, skip = LandLine, nrows=Ptc.list$PopN, col.names=PopLandInd_df_colnames )
    
    ## Landscape Indices - Summary
    ## Summary values of the Landscape Indices for each population
    LandIndSummary_colnames <- c("tot.edge","tot.core","avg.shp","avg.fractal","avg.perimeter")
    Ptc.list$LandIndSummary <- read.table( ptcFile, header=FALSE, skip = (LandLine+Ptc.list$PopN), nrows=1, col.names=LandIndSummary_colnames )

    ## Map averages
    ## This information is simply being stored as an unsorted list for the time being.
    MapAvgLine <- grep("Map averages",Ptc)
    Ptc.list$MapAvg <- Ptc[ (MapAvgLine+1):(MapAvgLine+Ptc.list$MapN) ]
    
  }
  
  return(Ptc.list)
} # End ptc.read function