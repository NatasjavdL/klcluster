#'@title simplify curve to given nr data points
#'
#'@description The function simplifies a polygonal curve P by an approximation
#' P'. Using a binary search on epsilon under the frechet distance, the desired
#' number of data points is returned.
#'
#'@param path This is a data.frame containing two columns with
#'information about the original trajectory/path. The columns are named x and y
#'and contain the x and y coordinates of the trajectory.
#'
#'@param size This is the number of data points the simplified path should have
#'
#'@param precision This is the error allowed on the number of data points. This
#'needs to be at least 1
#'
#'@export


bin_simplify <- function(path, size, precision) {

  if(size >= nrow(path)) {

    return(path)

  }

  # get upperbound
  upper <- frechet_distance(path, path[c(1,nrow(path)), ], 0.001)
  lower <- 0

  #init size
  path_size <- size + precision + 1

  #bin search on eps
  while(abs(path_size - size) > precision) {

    mid <- (lower + upper) / 2

    new_path <- simplify_curve(path, mid)

    path_size <- nrow(new_path)

    if(path_size > size) {

      lower <- mid

    } else {

      upper <- mid

    }

  }

  return(new_path)

}
