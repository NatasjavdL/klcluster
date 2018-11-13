# wrapper functions for Rcpp

#'@title calculate continuous frechet distance
#'
#'@description The function calculates the frechet distance between two paths.
#'The distance has a precision variable that indicates the precision of the
#'outcome.
#'
#'@param p1 This is a data.frame containing two columns with
#'information about the original trajectory/path. The columns are named x and y
#'and contain the x and y coordinates of the trajectory.
#'
#'@param p2 This is a data.frame containing two columns with
#'information about the original trajectory/path. The columns are named x and y
#'and contain the x and y coordinates of the trajectory.
#'
#'@param precision This is a number containing the precision with which to
#'calculate the frechet distance
#'
#'@return A number apx_eps which equals the frechet distance. Note that the
#'outcome is not exact. The true epsilon lies in the following range:
#' apx_eps - precision <= epsilon <= apx_eps
#'
#'@export

frechet_distance <- function(p1, p2, precision) {

  if(!is.matrix(p1)) {

    p1 <- df_to_mat(p1)
    p2 <- df_to_mat(p2)

  }

  return(get_frechet_c(p1, p2, precision))

}


#'@title compute DTW distance between two paths
#'
#'@param p1 This is a data.frame containing two columns with
#'information about the original trajectory/path. The columns are named x and y
#'and contain the x and y coordinates of the trajectory.
#'
#'@param p2 This is a data.frame containing two columns with
#'information about the original trajectory/path. The columns are named x and y
#'and contain the x and y coordinates of the trajectory.
#'
#'@param sqr_dist default = FALSE. If dist_measure = 3, then:
#'square_dist = TRUE ensures the dtw uses squared distances
#'square_dist = FALSE ensures the dtw uses normal distances
#'
#'@param normalize_dist default = TRUE. If dist_measure = 3, then:
#'norm_dist = TRUE ensures dtw distances returns a normalized distance
#' i.e. dtw_dist = dtw_dist/max(|p1|, |p2|)
#' if square_dist = TRUE then dtw_dist = sqrt(dtw_dist/max(|p1|, |p2|))
#'norm_dist = FALSE outputs actual dtw distance
#'
#'@export
#'
dtw_dist <- function(p1, p2, sqr_dist, normalize_dist) {

  # if paths not yet as matrix - transform
  if(!is.matrix(p1)) {

    p1 <- df_to_mat(p1)
    p2 <- df_to_mat(p2)

  }

  # get distance matrix dtw
  dtw_mat <- dtw(p1, p2, sqr_dist)

  # get final distance
  dist = dtw_mat[nrow(p1), nrow(p2)]

  # normalize if TRUE
  if(normalize_dist) {

    n_points = max(nrow(p1), nrow(p2))

    dist = dist/n_points

    if(sqr_dist) {

      dist = sqrt(dist)

    }

  }

  return(dist)

}


get_similarity <- function(method, p1, p2, param = 0.0001,
  sqr_dist = FALSE, norm_dist = TRUE) {
   # 1: cont frechet
   # 2: disc frechet
   # 3: dtw


  if(method == 1) {

    dist <- frechet_distance(p1, p2, param);

  } else if(method == 2) {

    dist <- discrete_frechet(p1, p2)

  } else {

    dist <- dtw_dist(p1, p2, sqr_dist, norm_dist)

  }


  return(dist)

}


