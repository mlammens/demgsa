#' @export metapopversion
#' @title Determine RAMAS Metapop version used to create a *.MP file
#' @description \code{metapopversion} is a function called from within the
#' \code{mp.read} function. It's single purpose is to read the `list` object
#' that contains all of the *.MP file values and determine which version
#' of RAMAS Metapop was used to construct the file.
#'
#' @param mpFile A list object of all of the parameters set in a *.MP file
#'
#'@return RAMAS Metapop version number as an integer
#'
metapopversion <- function(mpFile) {
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
