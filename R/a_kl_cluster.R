#'@title (k,l) clustering
#'
#'@description The function computes an number of centers of the given input
#'paths, where the centers are simplified to l data points. The number of
#'centers can be predefined when input k is given. However, the algorithm
#'can also determine the centers based on a given radius. When both k and
#'radius are specified, the algorithm returns a number of centers depending on
#'which condition is satisfied first
#'
#'@param paths This is a list of data.frames. Each data.frame represents a path.
#'Each data.frame contains two columns. The columns are named x and y
#'and contain the x and y coordinates of the trajectory
#'
#'@param radius The algorithm computes the number of centers such that each
#'path's distance to one of the centers is smaller than or equal to radius
#'
#'@param k The algorithm computs k centers
#'
#'@param l The number of data points to simplify each center to.
#'
#'@param prec_l This is a number containing the precision of the number of data
#'points to return.
#'
#'@param dist_measure either input 1, 2 or 3 correspondig to
#'"continuous frechet", "discrete frechet" or "dtw"
#'
#'@param prec_eps This is a number containing the precision with which to
#'calculate the frechet distance
#'
#'@param method Either "minmax" or "prob" for choosing centers
#'
#'@param square_dist default = FALSE. If dist_measure = 3, then:
#'square_dist = TRUE ensures the dtw uses squared distances
#'square_dist = FALSE ensures the dtw uses normal distances
#'
#'@param norm_dist default = TRUE. If dist_measure = 3, then:
#'norm_dist = TRUE ensures dtw distances returns a normalized distance
#' i.e. dtw_dist = dtw_dist/max(|p1|, |p2|)
#' if square_dist = TRUE then dtw_dist = sqrt(dtw_dist/max(|p1|, |p2|))
#'norm_dist = FALSE outputs actual dtw distance
#'
#'@return Clustering output containing clusters, centers and result
#'
#'@export

k_l_cluster <- function(paths, radius = NA, k = NA, l = 20, prec_l = 1,
  dist_measure = 2, prec_eps = 0.0001, method = "minmax", square_dist = FALSE,
  norm_dist = TRUE) {

  if(is.na(radius) && is.na(k)) {

    message("no radius or k given as input")

    return(NULL)

  } else if(is.na(radius) && !is.na(k)) {

    return(k_kl(paths, k, l, prec_l, dist_measure, prec_eps, method,
      square_dist, norm_dist))

  } else if(!is.na(radius) && is.na(k)) {

    return(radius_kl(paths, radius, l, prec_l, dist_measure, prec_eps, method,
      square_dist, norm_dist))

  } else {

    return(radius_k_kl(paths, radius, k, l, prec_l, dist_measure, prec_eps,
      method, square_dist, norm_dist))

  }

}

# kl clusterig based on radius

radius_kl <- function(paths, radius, l, prec_l, dist_measure, prec_eps, method,
  square_dist, norm_dist) {

  # Initialize algo:

  nr_paths <- length(paths)

  # create empty list T

  center_paths <- vector("list", nr_paths)

  # create empty list of ids
  path_ids_all <- vector("integer", nr_paths)

  # pick random path

  rand_path <- sample(nr_paths, 1)

  # add to path ids

  path_ids_all[1] <- rand_path

  # simplify path

  new_path <- bin_simplify(paths[[rand_path]], l, prec_l)

  # add simplified to set T
  center_paths[[1]] <- new_path

  # initialize distance table
  center_dist <- matrix(data = Inf, nr_paths, nr_paths)

  # fill first column distance matrix
  center_dist[, 1] <- sapply(1:nr_paths, function(x)
    get_similarity(dist_measure, new_path, paths[[x]], prec_eps, square_dist,
      norm_dist))

  # set radius
  rad <- max(center_dist[, 1])

  # check radius input
  if(radius >= rad) {

    message("radius too large, only one center returned")

  }

  if(radius < min(center_dist[, 1])) {

    message("radius very small, high probability of nr centers = nr paths")

  }

  # init counter
  nr_centers <- 1

  # start loop
  while(rad > radius && nr_centers < nr_paths) {

    nr_centers <- nr_centers + 1
    # find closest centers

    row_mins <- apply(center_dist, 1, min)

    # get new center

    if(method == "minmax") {

      path_id <- min_max(row_mins)

    } else if(method == "prob") {

      path_id <- sample_prob(row_mins)

    }

    # add to path id list

    path_ids_all[nr_centers] <- path_id

    # simplify

    center_paths[[nr_centers]] <- bin_simplify(paths[[path_id]], l, prec_l)

    # add new distances to table

    center_dist[, nr_centers] <- sapply(1:nr_paths, function(x)
      get_similarity(dist_measure, center_paths[[nr_centers]], paths[[x]],
        prec_eps, square_dist, norm_dist))

    # update radius
    row_mins <- apply(center_dist, 1, min)
    rad <- max(row_mins)

  }

  # trim of excess data
  path_ids_all <- path_ids_all[1:nr_centers]
  center_dist <- center_dist[, 1:nr_centers]
  center_paths <- center_paths[1:nr_centers]

  # create output
  row_mins <- apply(center_dist, 1, min)
  measures <- c("continuous frechet", "discrete frechet", "DTW")

  params <- list()
  params$radius <- NA
  params$k <- NA
  params$l <- l
  params$prec_l <- prec_l
  params$prec_eps <- ifelse(dist_measure == 1, yes = prec_eps, no = NA)
  params$dist_measure <- measures[dist_measure]

  results <- list()
  results$cost <- max(row_mins)
  results$minsum <- sum(row_mins)
  results$center_ids <- path_ids_all
  results$distances <- center_dist

  # set names center dist
  colnames(center_dist) <- paste0("p", path_ids_all)

  # assign all info to list structure
  out <- list()

  out$centers <- center_paths
  out$clusters <- get_clusters(center_dist)
  out$result <- results
  out$parameters <- params
  # return center_paths
  return(out)


}

# kl clustering based on radius and k

radius_k_kl <- function(paths, radius, k, l, prec_l, dist_measure, prec_eps,
  method, square_dist, norm_dist) {

  # Initialize algo:

  nr_paths <- length(paths)

  # create empty list T

  center_paths <- vector("list", k)

  # create empty list of ids
  path_ids_all <- vector("integer", k)

  # pick random path

  rand_path <- sample(nr_paths, 1)

  # add to path ids

  path_ids_all[1] <- rand_path

  # simplify path

  new_path <- bin_simplify(paths[[rand_path]], l, prec_l)

  # add simplified to set T
  center_paths[[1]] <- new_path

  # initialize distance table
  center_dist <- matrix(data = Inf, nr_paths, k)

  center_dist[, 1] <- sapply(1:nr_paths, function(x)
    get_similarity(dist_measure, new_path, paths[[x]], prec_eps, square_dist,
      norm_dist))

  # init radius and i (k)
  i <- 1
  rad <- max(center_dist[, 1])

  # check radius input
  if(radius >= rad) {

      message("radius and k too large, only one center returned")

  }

  if(radius < min(center_dist[, 1])  || k >= nr_paths) {

    if(k < nr_paths) {

      message("radius very small, algorithm will rely on k")

    } else {

      message("radius very small and k too large,
        high probability of nr centers = nr paths")

    }

  }

  # start loop
  if(k > 1) {

    while(i < k && rad > radius) {

      # update i
      i <- i + 1

      # find closest centers

      row_mins <- apply(center_dist, 1, min)

      # get new center

      if(method == "minmax") {

        path_id <- min_max(row_mins)

      } else if(method == "prob") {

        path_id <- sample_prob(row_mins)

      }

      # add to path id list

      path_ids_all[i] <- path_id

      # simplify

      center_paths[[i]] <- bin_simplify(paths[[path_id]], l, prec_l)

      # add new distances to table

      center_dist[, i] <- sapply(1:nr_paths, function(x)
        get_similarity(dist_measure, center_paths[[i]], paths[[x]], prec_eps,
          square_dist, norm_dist))



      # update radius
      row_mins <- apply(center_dist, 1, min)
      rad <- max(row_mins)

    }

    # trim of excess data

    if(i < k) {

      path_ids_all <- path_ids_all[1:i]
      center_dist <- center_dist[, 1:i]
      center_paths <- center_paths[1:i]

    }

  }


  # create output
  row_mins <- apply(center_dist, 1, min)
  measures <- c("continuous frechet", "discrete frechet", "DTW")

  params <- list()
  params$radius <- radius
  params$k <- k
  params$l <- l
  params$prec_l <- prec_l
  params$prec_eps <- ifelse(dist_measure == 1, yes = prec_eps, no = NA)
  params$dist_measure <- measures[dist_measure]

  results <- list()
  results$cost <- max(row_mins)
  results$minsum <- sum(row_mins)
  results$center_ids <- path_ids_all
  results$distances <- center_dist

  # set names center dist
  colnames(center_dist) <- paste0("p", path_ids_all)

  # assign all info to list structure
  out <- list()

  out$centers <- center_paths
  out$clusters <- get_clusters(center_dist)
  out$result <- results
  out$parameters <- params


  # return center_paths

  return(out)

}


# kl clustering based on k

k_kl <- function(paths, k, l, prec_l, dist_measure, prec_eps, method,
  square_dist, norm_dist) {

  # Initialize algo:

  nr_paths <- length(paths)

  # check value k
  if(k >= nr_paths) {

    message("k too large for data")
    return(NULL)

  }

  # create empty list T

  center_paths <- vector("list", k)

  # create empty list of ids
  path_ids_all <- vector("integer", k)

  # pick random path

  rand_path <- sample(nr_paths, 1)

  # add to path ids

  path_ids_all[1] <- rand_path

  # simplify path

  new_path <- bin_simplify(paths[[rand_path]], l, prec_l)

  # add simplified to set T
  center_paths[[1]] <- new_path

  # initialize distance table
  center_dist <- matrix(data = Inf, nr_paths, k)

  center_dist[, 1] <- sapply(1:nr_paths, function(x)
    get_similarity(dist_measure, new_path, paths[[x]], prec_eps, square_dist,
      norm_dist))

  # start loop
  if(k > 1) {

    for(i in 2:k) {

      pids <- setdiff(1:length(paths), path_ids_all)

      # find closest centers

      row_mins <- apply(center_dist[pids, ], 1, min)

      # get new center

      if(method == "minmax") {

        path_id <- pids[min_max(row_mins)]

      } else if(method == "prob") {

        path_id <- pids[sample_prob(row_mins)]

      }

      # add to path id list

      path_ids_all[i] <- path_id

      # simplify

      center_paths[[i]] <- bin_simplify(paths[[path_id]], l, prec_l)

      # add new distances to table

      center_dist[, i] <- sapply(1:nr_paths, function(x)
        get_similarity(dist_measure, center_paths[[i]], paths[[x]], prec_eps,
          square_dist, norm_dist))


    }

  }

  # create output
  row_mins <- apply(center_dist, 1, min)
  measures <- c("continuous frechet", "discrete frechet", "DTW")

  params <- list()
  params$radius <- NA
  params$k <- NA
  params$l <- l
  params$prec_l <- prec_l
  params$prec_eps <- ifelse(dist_measure == 1, yes = prec_eps, no = NA)
  params$dist_measure <- measures[dist_measure]

  results <- list()
  results$cost <- max(row_mins)
  results$minsum <- sum(row_mins)
  results$center_ids <- path_ids_all
  results$distances <- center_dist

  # set names center dist
  colnames(center_dist) <- paste0("p", path_ids_all)

  # assign all info to list structure
  out <- list()

  out$centers <- center_paths
  out$clusters <- get_clusters(center_dist)
  out$result <- results
  out$parameters <- params
  # return center_paths

  return(out)

}


