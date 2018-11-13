#'@title kl-median
#'
#'@description The function computes k centers of the given input paths and
#'initial centers, where the centers are simplified to l data points.
#'
#'@param paths This is a list of data.frames. Each data.frame represents a path.
#'Each data.frame contains two columns. The columns are named x and y
#'and contain the x and y coordinates of the trajectory
#'
#'@param clust_out the output of one of the clustering algorithms
#'
#'@param max_iter integer with maximum number of iterations to perform
#'
#'@param dist_measure either input 1, 2 or 3 correspondig to
#'"continuous frechet", "discrete frechet" or "dtw"
#'
#'@param clusts_to_use (optional)vector of integers with indexes of clusters
#'to use in the algorithm#'
#'
#'@param prec_eps This is a number containing the precision with which to
#'calculate the frechet distance
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
#'@param update if TRUE then all output will be updated
#'
#'@return Clustering output containing clusters, centers and result
#'
#'@export

kl_median <- function(paths, clust_out, max_iter = 1, dist_measure = 2,
  clusts_to_use = NA, prec_eps = 0.0001, square_dist = FALSE,
  norm_dist = TRUE, update = TRUE) {

  if(is.na(clusts_to_use)) {

    clusts_to_use <- 1:length(clust_out$clusters)

  }

  if(length(clusts_to_use) == 1) {

    update <- FALSE
    max_iter <- 1

  }

  improvement <- TRUE

  new_out <- clust_out

  for(i in 1:max_iter) {

    for(j in 1:length(clusts_to_use)) {

      clust <- clusts_to_use[j]
      center <- new_out$centers[[clust]]
      center_id <- new_out$result$center_ids[clust]
      path_ids <- new_out$clusters[[clust]]

      new_centers <- kl_median_cluster(paths, center, path_ids,
        center_id, l = clust_out$parameters$l, clust_out$parameters$prec_l,
        dist_measure, prec_eps, square_dist, norm_dist)

      new_out$centers[[clust]] <- new_centers[[1]]
      new_out$result$center_ids[clust] <- new_centers[[2]]


    }

    # get distance matrix
    distm <- data.frame("c0" = 1:length(paths))

    for(i in 1:length(new_out$centers)) {

      distm[paste0("c", i)] <- sapply(1:length(paths), function(x)
        get_similarity(dist_measure, new_out$centers[[i]], paths[[x]], prec_eps))

    }

    distm["c0"] = NULL

    row_mins <- apply(distm, 1, min)

    results <- list()
    results$cost <- max(row_mins)
    results$minsum <- sum(row_mins)
    results$center_ids <- new_out$result$center_ids
    results$distances <- df_to_mat(distm)


    if(update) { new_out$clusters <- get_clusters(distm) }
    new_out$result <- results


    if(new_out$result$minsum < clust_out$result$minsum) {

      message(paste0("Improvement minsum from ", clust_out$result$minsum,
        " to ", new_out$result$minsum))

      clust_out <- new_out

    } else {

      message("No further improvements made")

      return(new_out)

    }

  }

  return(new_out)

}


kl_median_cluster <- function(paths, center, path_ids, center_id, l, prec_l,
  dist_measure, prec_eps, square_dist, norm_dist) {

  # Initialize algo:

  nr_paths <- length(path_ids)

  paths_to_use <- paths[path_ids]
  paths_min_center <- setdiff(path_ids, center_id)

  curr_dist <- sum(sapply(1:nr_paths, function(x)
    get_similarity(dist_measure, center, paths_to_use[[x]], prec_eps,
      square_dist, norm_dist)))

  for(i in paths_min_center) {

    simple_path <- bin_simplify(paths[[i]], l, prec_l)

    new_dist <- sum(sapply(1:nr_paths, function(x)
      get_similarity(dist_measure, simple_path, paths_to_use[[x]], prec_eps,
        square_dist, norm_dist)))

    if(new_dist < curr_dist) {

      curr_dist <- new_dist
      center <- simple_path
      center_id <- i

    }

  }

  return(list(center, center_id))

}
