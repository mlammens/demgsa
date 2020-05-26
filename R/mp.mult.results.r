#' @export mp.mult.results
#' @title Extract results for a list of *.MP files
#' @description
#' Call \code{mp.results} for a list of files and export results as a csv file
#'
#' @param mp.file.list A list (character vector) of *.MP files to extract
#' results from
#' @param out.csv A name of a csv file to write results to
#' @param spatial (TRUE/FALSE) Do the *.MP files contain spatial information
#' to be extracted
#'
#' @seealso \code{\link{mp.results}} for a full list of descriptions for each
#' of the columns in the resulting CSV file.

mp.mult.results <- function( mp.file.list, out.csv, spatial=FALSE ) {
  # Call mp.results for a list of files and export results as a csv file
  #
  # Author: Matthew Aiello-Lammens
  # Created: 23 January 2011
  #
  # Args:
  #  mp.file.list: a text file of the list of *.mp files for which results should be extracted
  #
  # Returns:
  #  No direct returns. Creates a CSV file of input and output metrics for multiple
  #  *.mp files
  #
  #
  ###################################################################################################

  # mp.list <- readLines( mp.file.list )
  mp.list <- mp.file.list
  mp.file.cnt <- length( mp.list )

  mp.mult.res <- vector()
  for ( mp in 1:mp.file.cnt ) {
    mp.res <- mp.results( mp.list[mp], spatial=spatial )
    mp.mult.res <- rbind( mp.mult.res, mp.res )
  }
  ## Rename first column
  names( mp.mult.res )[1] <- "mpFile"
  # Write whole table after all files have been combined.
  write.csv( mp.mult.res, out.csv )
}
