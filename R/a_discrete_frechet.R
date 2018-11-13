#'@title compute discrete frechet distance
#'
#'@description from Thomas Eiter and Heikki Mannila
#'
#'@param p1 data.frame containing x and y coordinates of trajectory
#'
#'@param p2 data.frame containing x and y coordinates of trajectory
#'
#'@references \url{http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf}
#'
#'@export
#'

discrete_frechet <- function(p1, p2) {

  if(!is.matrix(p1)) {

    p1 <- df_to_mat(p1)
    p2 <- df_to_mat(p2)

  }

  dist_mat <- matrix(-1, nrow(p1), nrow(p2))

  dist <- calc_ca(dist_mat, nrow(p1) - 1, nrow(p2) - 1, p1, p2)

  return(dist)

}
