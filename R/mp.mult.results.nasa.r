mp.mult.results.nasa <- function( mp.file.list, out.csv, spatial=TRUE ) {
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

# Get data about matching Species and Generic Demographic Model
load(paste(sens.base.dir,'spec.mod.RData',sep=''))
  
mp.list <- readLines( mp.file.list )
mp.file.cnt <- length( mp.list )

# NASA Specific Column Names
mp.res.colnames <- c("file","spec","model","clim.chg","mp.run","match.stat","metapop.initab","Gen.Abar",
                     "Gen.Mu1","Gen.Tgen","EigenVal","dd.type","Rmax","GrowthRt","fec.stdev.avg","surv.stdev.avg",
                     "stdev.avg","fec.cv.avg","surv.cv.avg","cv.avg","N.maxVmin.10","N.plusVmin.SD","N.CV.10",
                     "St.First.Rep","St.First.Rep.Name","ext.risk","threshold","prob.thresh.maxt","med.time.cross",
                     "prob.50","prob.250","prob.1000","prob.Thr.50","exp.min.n","sderr.ema","n.mean","n.stdev",
                     "metapop.chng","occ.maxt","occ.stdev.maxt","quant.05","quant.25","quant.50","quant.75",
                     "quant.95","harv.avg","harv.stdev","harv.min","harv.max","use.disp.dist.func",
                     "disp.dist.func.a",
                     "avg.disp.dist.b","max.disp.dist.Dmax","mean.t0.disp.rate","loquart.t0.disp.rate",
                     "hiquart.t0.disp.rate","mean.t0.Ngt0.disp.rate","loquart.t0.Ngt0.disp.rate",
                     "hiquart.t0.Ngt0.disp.rate","mean.t0.nearest.disp.rate","loquart.t0.nearest.disp.rate",
                     "hiquart.t0.nearest.disp.rate","mean.t0.Ngt0.nearest.disp.rate",
                     "loquart.t0.Ngt0.nearest.disp.rate","hiquart.t0.Ngt0.nearest.disp.rate","use.corr.dist.func",
                     "avg.corr.dist.b",
                     "metapop.abundance.1","avg.dist.occupied.pops.1",
                     "avg.dist.near.neighbor.1","low.quart.dist.near.neighbor.1",
                     "upper.quart.dist.near.neighbor.1","max.dist.near.neighbor.1","frag.index.hra.1",
                     "patch.n.1","tot.patch.area.1",
                     "tot.core.area.1","tot.edge.length.1","lg.patch.area.1","lg.patch.core.1",
                     "overall.edge.area.ratio.1","overall.shape.ind.1","overall.frac.dim.1","avg.patch.area.1",
                     "avg.core.area.1","avg.edge.length.1","avg.shape.ind.1","avg.frac.dim.1",
                     "metapop.abundance.10","avg.dist.occupied.pops.10",
                     "avg.dist.near.neighbor.10","low.quart.dist.near.neighbor.10",
                     "upper.quart.dist.near.neighbor.10","max.dist.near.neighbor.10","frag.index.hra.10",
                     "patch.n.10","tot.patch.area.10",
                     "tot.core.area.10","tot.edge.length.10","lg.patch.area.10","lg.patch.core.10",
                     "overall.edge.area.ratio.10","overall.shape.ind.10","overall.frac.dim.10","avg.patch.area.10",
                     "avg.core.area.10","avg.edge.length.10","avg.shape.ind.10","avg.frac.dim.10",
                     "metapop.abundance.11","avg.dist.occupied.pops.11",
                     "avg.dist.near.neighbor.11","low.quart.dist.near.neighbor.11",
                     "upper.quart.dist.near.neighbor.11","max.dist.near.neighbor.11","frag.index.hra.11",
                     "patch.n.11","tot.patch.area.11",
                     "tot.core.area.11","tot.edge.length.11","lg.patch.area.11","lg.patch.core.11",
                     "overall.edge.area.ratio.11","overall.shape.ind.11","overall.frac.dim.11","avg.patch.area.11",
                     "avg.core.area.11","avg.edge.length.11","avg.shape.ind.11","avg.frac.dim.11",
                     "metapop.abundance.21","avg.dist.occupied.pops.21",
                     "avg.dist.near.neighbor.21","low.quart.dist.near.neighbor.21","upper.quart.dist.near.neighbor.21",
                     "max.dist.near.neighbor.21","frag.index.hra.21",
                     "patch.n.21","tot.patch.area.21","tot.core.area.21",
                     "tot.edge.length.21","lg.patch.area.21","lg.patch.core.21","overall.edge.area.ratio.21",
                     "overall.shape.ind.21","overall.frac.dim.21","avg.patch.area.21","avg.core.area.21",
                     "avg.edge.length.21","avg.shape.ind.21","avg.frac.dim.21",
                     "metapop.abundance.31","avg.dist.occupied.pops.31",
                     "avg.dist.near.neighbor.31","low.quart.dist.near.neighbor.31","upper.quart.dist.near.neighbor.31",
                     "max.dist.near.neighbor.31","frag.index.hra.31",
                     "patch.n.31","tot.patch.area.31","tot.core.area.31",
                     "tot.edge.length.31","lg.patch.area.31","lg.patch.core.31","overall.edge.area.ratio.31",
                     "overall.shape.ind.31","overall.frac.dim.31","avg.patch.area.31","avg.core.area.31",
                     "avg.edge.length.31","avg.shape.ind.31","avg.frac.dim.31",
                     "metapop.abundance.111","avg.dist.occupied.pops.111",
                     "avg.dist.near.neighbor.111","low.quart.dist.near.neighbor.111","upper.quart.dist.near.neighbor.111",
                     "max.dist.near.neighbor.111","frag.index.hra.111",
                     "patch.n.111","tot.patch.area.111","tot.core.area.111",
                     "tot.edge.length.111","lg.patch.area.111","lg.patch.core.111","overall.edge.area.ratio.111",
                     "overall.shape.ind.111","overall.frac.dim.111","avg.patch.area.111","avg.core.area.111",
                     "avg.edge.length.111","avg.shape.ind.111","avg.frac.dim.111")

write.table(t(mp.res.colnames),file=out.csv,col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)

mp.mult.res <- vector()
for ( mp in 1:mp.file.cnt ) {
  print(paste('******** Extracting Data from:',mp.list[mp] ,'********'))
  # NASA Specific data extraction - assumes a very specific file path format
  run.info <- unlist(strsplit( mp.list[mp], split='/'))
  rl <- length(run.info)
  mp.run <- run.info[ rl ]
  clim.ch <- run.info[ rl-1 ]
  mod <- run.info[ rl-2 ]
  spec <- run.info[ rl-3 ]
  # Determine if species and generic model are matching or not
  # First determine what the matching generic model should be for this species
###browser()
  spec.match <- spec.mod$Model[which(spec.mod$Species == spec)]
  # Test for match between species and generic model used
  match.stat <- ifelse( spec.match==mod, 'Match','Mis-match')
  run.info.df <- data.frame( spec=spec, model=mod, clim.chg=clim.ch, mp.run=mp.run, match.stat=match.stat )
  ## Create arguments for mp.results with spatial data, ptc files, and habdyn files.
  # Get ptcFiles and hdhFile
  
  # Swap out Population for Spatial in path name
  spatial.path <- sub('Population','Spatial',mp.list[mp])
  # Remove the generic model name
  spatial.path <- sub( mod, '', spatial.path)
  spatial.path <- sub( '//','/',spatial.path)
  
  if ( clim.ch == 'NoCC' ){
    ptc.name.1990 <- paste( spec, 'LO1990.ptc', sep='')
    ptc.name.full.1990 <- sub( mp.run, ptc.name.1990, spatial.path)
    ## ptcFiles <- c( ptc.name.full.1990, ptc.name.full.1990 ) #Worked previously
    ptcFiles <- c( ptc.name.full.1990, ptc.name.full.1990, ptc.name.full.1990, ptc.name.full.1990, ptc.name.full.1990, ptc.name.full.1990 )
    UseHabDyn <- FALSE
    hdhFile <- 'no file'
  } else {
    ptc.name.1990 <- paste( spec, 'LO1990.ptc', sep='')
    ptc.name.full.1990 <- sub( mp.run, ptc.name.1990, spatial.path)
    ptc.name.2010 <- paste( spec, 'LO2010.ptc', sep='')
    ptc.name.full.2010 <- sub( mp.run, ptc.name.2010, spatial.path)
    ## Below two ptc.name.full.2XXX added on 31 July 2012
    ptc.name.1999 <- paste( spec, 'LO1999.ptc', sep='')
    ptc.name.full.1999 <- sub( mp.run, ptc.name.1999, spatial.path)
    ptc.name.2000 <- paste( spec, 'LO2000.ptc', sep='')
    ptc.name.full.2000 <- sub( mp.run, ptc.name.2000, spatial.path)
    ptc.name.2020 <- paste( spec, 'LO2020.ptc', sep='')
    ptc.name.full.2020 <- sub( mp.run, ptc.name.2020, spatial.path)
    ptc.name.2100 <- paste( spec, 'LO2100.ptc', sep='')
    ptc.name.full.2100 <- sub( mp.run, ptc.name.2100, spatial.path)
    ##
    ptcFiles <- c(ptc.name.full.1990, ptc.name.full.1999, ptc.name.full.2000, ptc.name.full.2010, ptc.name.full.2020, ptc.name.full.2100 )
    UseHabDyn <- TRUE
    # Use some regexs to name the habdyn history file (hdhFile)
    hdhFileName <- paste( spec, 'LO-hist.TXT', sep='' )
    hdhFile <- sub( mp.run, hdhFileName, spatial.path )
  }
  
  # Call mp.results on individual *.mp file
  mp.res <- mp.results( mp.list[mp], spatial=spatial, mac=FALSE, ptc=TRUE, ptcFiles=ptcFiles,
                        ptcFileIter=c(1,10,11,21,31,111), habdyn=UseHabDyn, hdhFile=hdhFile)
  mp.res <- cbind( run.info.df, mp.res )
  write.table( mp.res, file=out.csv, append=TRUE, col.names=FALSE,sep=",")
  ###mp.mult.res <- rbind( mp.mult.res, mp.res )
}

###write.csv( mp.mult.res, out.csv )
}
