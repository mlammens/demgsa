sensitivity.setup <- function(sens.base.dir) {
# RAMAS sensitivity analysis program Setup
# Author: Matthew Aiello-Lammens
#
# Source additional R scripts 
# 
# Arguments:
#	sens.base.dir: The base directory where all of the sensitivity scripts reside.  
# 	Some examples of sens.base.dir are:
#		sens.base.dir='/Users/mlammens/Documents/My Dropbox/RA-Sensitivity/' # Sensitivity base directory for MA-L
#		sens.base.dir='/Users/mlammens/Dropbox/RA-Sensitivity/'
#		sens.base.dir='/Users/Matthew/Documents/My Dropbox/RA-Sensitivity/'
#
# Returns:
#	No returns frome this script
#
source( paste(sens.base.dir,'mp.read.r',sep="") )
source( paste(sens.base.dir,'mp.read.results.r',sep="") )
source( paste(sens.base.dir,'mp.write.r',sep="") )
source( paste(sens.base.dir,'ptc.read.r',sep="") )
source( paste(sens.base.dir,'hdhist.read.r',sep="") )
source( paste(sens.base.dir,'fill.matrix.df.r',sep="") )
source( paste(sens.base.dir,'metapopversion.r',sep="") )
source( paste(sens.base.dir,'sensitivity.r',sep="") )
source( paste(sens.base.dir,'mp.results.r',sep="") )
source( paste(sens.base.dir,'mp.mult.results.r',sep="") )
# Make sens.base.dir a global variable
sens.base.dir <<- sens.base.dir
# Load a char vector called 'guide', which is a guide to the mp.read list structure.  Here it is
# used to inform the user which parameters differ between *.mp files
load( paste(sens.base.dir,'mp.guide.Rdata',sep=""), .GlobalEnv )
}
# Load necessary packages
library('fields')
library('lhs')
library('popbio')
