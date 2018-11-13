#'@title simplify curves to given nr data points
#'
#'@description The function simplifies all polygonal curve P, stored in a list
#' by an approximation P'. Using a binary search on epsilon under the frechet
#' distance, the desired number of data points is returned for each path.
#'
#'@param paths This is a list of data.frames. Each data.frame contains two
#'columns with information about the original trajectory/path. The columns are
#'named x and y and contain the x and y coordinates of the trajectory.
#'
#'@param size This is the number of data points the simplified path should have
#'
#'@param precision This is the error allowed on the number of data points. This
#'needs to be at least 1
#'
#'@export
#'

simplify_all <- function(paths, size, precision) {

  # initialize simple paths
  simple_paths <- list()

  for(i in 1:length(paths)) {

    simple_paths[[i]] <- bin_simplify(paths[[i]], size, precision)

  }

  return(simple_paths)

}
