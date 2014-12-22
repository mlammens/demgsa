gen.time <- function( mp.list ) {
	gen.list <- vector()
	for ( l in mp.list ) {
		gen.call <- paste( '../GenTime ', '"', l , '"' )
		print(gen.call) ### DEBUG LINE
		gen <- system( gen.call, intern = TRUE )
		gen <- strsplit( gen, " ")
		gen <- unlist( lapply( gen, as.numeric ) )
		gen.list <- rbind(gen.list, gen)
	}
	cols <- c("Abar","Mu1","Tgen","domEigenvalue","errorcode")
	gen.df <- data.frame(gen.list, row.names=mp.list)
	names(gen.df) <- cols
	return(gen.df)
}