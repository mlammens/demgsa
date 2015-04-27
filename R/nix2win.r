#' @export nix2win
#' @title Convert MP files from *nix to Windows format
#' @description
#' \code{nix2win} is a wrapper function for a bash script to convert from
#' *nix formated text files to Windows formated text files. These files have
#' different end-of-line characters. RAMAS GIS only reads Windows formated
#' files.
#'
#' @param mp.new.file Name of the *.mp file to convert from *nix to Windows
#' format
#'
#' @return
#' There are no returns to the R environment. However, the selected *.mp file
#' will be in Windows format after the function has been executed.


# A simple script to convert a unix file to Windows format

nix2win <- function( mp.new.file ){
  ## Define the *nix shell script needed to convert a file from *nix to windows
  ## format
  nix2win.shell <- paste0( "#!/bin/bash\n",
                           "f=",mp.new.file,"\n",
                           #cat(paste0( "f_nix=", '"', "${f}_nix", '"', "\n" )),
                           'f_nix="${f}_nix"\n',
                           "mv $f $f_nix\n",
                           #awk 'sub("$", "\r")' $f_nix > $f
                           'awk \'sub("$", "\\r")\' $f_nix > $f\n',
                           "rm $f_nix" )

  system( nix2win.shell )
}
