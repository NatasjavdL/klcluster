# script for .onLoad() and .onAttach()

#' @useDynLib klcluster
#' @importFrom Rcpp sourceCpp
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("klcluster", libpath)
}
