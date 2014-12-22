mp.mult.results <- function( mp.file.list, out.csv, spatial=FALSE, mac=FALSE ) {
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

mp.list <- readLines( mp.file.list )
mp.file.cnt <- length( mp.list )

mp.mult.res <- vector()
for ( mp in 1:mp.file.cnt ) {
  mp.res <- mp.results( mp.list[mp], spatial=spatial, mac=mac )
  mp.mult.res <- rbind( mp.mult.res, mp.res )
  ###write.table( mp.res, file=out.csv, append=TRUE, col.names=FALSE,sep=",")
}
# Write whole table after all files have been combined. 
write.csv( mp.mult.res, out.csv )
}
