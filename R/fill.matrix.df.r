fill.matrix.df <- function(PopDF, Func, Type) {
	# Creates a Dispersal or Correlation Matrix based on the dispersal-(correlation-)distance function
	#
	# Args:
	#	PopList: Population information from the mp.file structure; specifically,
	#	the information in the mp.file$PopList element
	#
	#	Func: Dispersal- or Correlation-Distance function parameters.  Element mp.file$DispDistFunc
	#	or mp.file$CorrDistFunc in the mp.file 'list' structure.
	#
	#	Type: The type of matrix to fill (either 'disp' for dispersal or 'corr' for correlation)
	#
	# Returns:
	#	DispMatr: A dispersal matrix
	#
	# Author: Matthew Aiello-Lammens
	# Date: 11 February 2011
	###############################################################################################
	if ( Type == 'disp' | Type == 'corr' ) {
		PopNumber <- nrow(PopDF)
		pop.coords <- cbind(PopDF$X_coord,PopDF$Y_coord)
		DistMatr <- rdist(pop.coords)
		if ( Type == 'disp' ) {
			# Func variable is a list of dispersal-distance function parameters. The function is
			# Mij = a exp(-Dij^c/b). Func = c(a, b, c, Dmax).  Dmax is the maximum dispersal rate.
			DispMatr <- Func[1] * exp( (-DistMatr^Func[3]) / Func[2] )
			diag(DispMatr) <- 0 # Replace diagnols with 0's
			# If the distance is greater than Dmax, make dispersal equal to 0
			DispMatr[ which( DistMatr > Func[4] ) ] <- 0
			# If the dispersal value is less than 0.0000001, set equal to 0
			DispMatr[ which( DispMatr < 0.0000001 ) ] <- 0
			return(DispMatr)
		} else if ( Type == 'corr' ) {
			# Func variable is a list of corr-distance function parameters. The function is
			# Cij = a exp(-Dij^c/b). Func = c(a, b, c).
			CorrMatr <- Func[1] * exp( (-DistMatr^Func[3]) / Func[2] )
			# Make the matrix only lower triangular
			CorrMatr <- lower.tri(CorrMatr) * CorrMatr
			# If the correlation is less than 0.0000001, set equal to 0
			CorrMatr[ which( CorrMatr < 0.0000001 ) ] <- 0
			diag(CorrMatr) <- 1 # Replace diagnols with 1's
			return(CorrMatr)
		}
	} else {
		stop("Incorrect 'Type' called in fill.matrix function.  Type must be either 'disp' or 'corr'")
	}
} # End fill.matrix function
