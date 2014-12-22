metapopversion <- function(mpFile) {
	# Reads a 'list' variable that is composed of all lines of a *.mp file, as used
	# in the RAMAS Metapop program, and returns the Metapop version number as an integer.
	#
	# Args:
	#	mpFile: a 'list' variable that is composed of all of the lines of a *.mp file
	#
	# Returns:
	#	RAMAS Metapop version number, as an integer
	#
	# Author: Matthew Aiello-Lammens
	# Date: 2 February 2011
	
	verLine <- mpFile[1]
	verLine <- unlist(strsplit(verLine,' '))
	if (verLine[1] == "Metapopulation") {
		version <- verLine[4] # Version number should be the fourth element
		version <- gsub("[^0-9]","",version) # Remove non-numeric values
		return(version)
	} else {
		stop("Error: *.mp File is not a RAMAS Metapop File")
	}
}
