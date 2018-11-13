#'@title single linkage clustering
#'
#'@description This algorithm builds a hierarchy of clusterings and returns
#'a number of centers given either radius or k
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
#'@export
#'

single_link <- function(paths, radius = NA, k = NA, l = 20, prec_l = 1,
  dist_measure = 2, prec_eps = 0.0001, square_dist = FALSE, norm_dist = TRUE) {

  if(is.na(radius) && is.na(k)) {

    message("no radius or k given as input")

    return(NULL)

  }

  if(!is.list(paths)) {

    message("paths are in incorrect format")

    return(NULL)

  }

  # Initialize algo:

  nr_paths <- length(paths)

  # create distance matrix to define clusters

  distm <- matrix(NA_integer_, nr_paths, nr_paths)

  for(i in 1:(nr_paths - 1)) {

    j <- i + 1

    tmp <- sapply(j:nr_paths, function(x) get_similarity(dist_measure,
      paths[[i]], paths[[x]], square_dist, norm_dist))

    distm[j:nr_paths, i] <- tmp
    distm[i, j:nr_paths] <- tmp

  }

  # set diag to zero temporarily
  diag(distm) <- 0

  # keep copy original for finding centers & radius
  distm_copy <- distm

  #### define method to use

  if(is.na(radius) && is.na(k)) {

    message("no radius or k given as input")

    return(NULL)

  } else if(is.na(radius) && !is.na(k)) {

    clusterm <- k_single(distm, nr_paths, k)

  } else if(!is.na(radius) && is.na(k)) {

    if(max(distm) < radius) {

      message("radius too large")
      return(NULL)

    }

    clusterm <- radius_single(distm, distm_copy, nr_paths, radius)

  } else {

    clusterm <- radius_k_single(distm, distm_copy, nr_paths, radius, k)

  }

  # set cluster numbering for single trajectories
  if(is.na(sum(clusterm[, 2]))) {

    for(i in 1:nr_paths) {

      if(is.na(clusterm[i, 2])) {

        clusterm[i, 2] <- max(clusterm[, 2], na.rm = TRUE) + 1

      }

    }

  }


  # adjust cluster numbering if not sequential
  clust_vec <- clusterm[, 2]
  clust_ids <- unique(clusterm[, 2])
  max_clust <- max(clust_ids)

  if(max_clust > length(clust_ids)) {

    tmp <- rep(NA_integer_, max_clust)
    tmp[clust_ids] <- 1:length(clust_ids)
    clust_vec <- mapply(function(x) x = tmp[x], x = clusterm[, 2])

  }

  # create list of clusters

  clust_list <- list()

  for(i in 1:max(clust_vec)) {

    clust_list[[i]] <- which(clust_vec == i)

  }

  # find centers

  # initialize
  path_ids_all <- vector(mode = "integer", length(clust_list))
  center_dist <- matrix(data = Inf, nr_paths, length(clust_list))
  center_paths <- list()

  for(i in 1:length(clust_list)) {

    # get cluster
    cluster <- clust_list[[i]]

    # check if more than 2 paths
    if(length(cluster) < 3) {

      path_id <- cluster[1]


    } else {

      dm <- distm_copy[cluster, cluster]

      diag(dm) <- 0

      # use minimum of sum of distaces to determine center
      center_id <- which(colSums(dm) == min(colSums(dm)))[1]

      path_id <- cluster[center_id]

    }

    path_ids_all[i] <- path_id
    center <- paths[[path_id]]
    center <- bin_simplify(center, l, prec_l)
    center_paths[[i]] <- center

    center_dist[, i] <- sapply(1:nr_paths, function(x)
      get_similarity(dist_measure, center, paths[[x]], prec_eps, square_dist,
        norm_dist))

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
  out$clusters <- clust_list
  out$result <- results
  out$parameters <- params


  # return center_paths

  return(out)

}


#### single linkage radius based ####

radius_single <- function(distm, distm_copy, nr_paths, radius) {

  # get rid of distances to paths itself - inf

  diag(distm) <- Inf

  # define cluster matrix - col1 = paths, col2 = cluster id

  clusterm <- matrix(0, nr_paths, 2)

  clusterm[, 1] <- 1:nr_paths
  clusterm[, 2] <- NA

  # find clusters based on radius

  #initialize cluster indexing and radius
  nr_clusts <- nr_paths
  clust_counter <- 1
  rad <- 0

  while(rad < radius) {

    min_index <- which(distm == min(distm), arr.ind = TRUE)

    p1 <- min_index[1,1]
    p2 <- min_index[1,2]


    # both not assigned to cluster yet
    if(is.na(clusterm[p1, 2]) && is.na(clusterm[p2, 2])) {

      clusterm[p1, 2] <- clust_counter
      clusterm[p2, 2] <- clust_counter

      distm[p1, p2] <- Inf
      distm[p2, p1] <- Inf

      # update counter nr clusters
      clust_counter <- clust_counter + 1

      # update nr clusters
      nr_clusts <- nr_clusts - 1


    } else if(!is.na(clusterm[p1, 2]) && !is.na(clusterm[p2,2])) {

      # if both already in cluster

      id1 <- clusterm[p1, 2]
      id2 <- clusterm[p2, 2]

      # merge clusters
      clusterm[which(clusterm[,2] == id2), 2] <- id1

      clust <- which(clusterm[,2] == id1)

      # set distm to infinity to rule out inner cluster distances
      distm[clust, clust] <- Inf

      # update nr of clusters
      nr_clusts <- nr_clusts - 1

    } else if(is.na(clusterm[p1, 2])) {

      # if p2 already in cluster

      clusterm[p1, 2] <- clusterm[p2, 2]

      clust <- which(clusterm[,2] == clusterm[p2, 2])

      # set distm to infinity to rule out inner cluster distances
      distm[clust, clust] <- Inf

      # update nr clusters
      nr_clusts <- nr_clusts - 1

    } else {

      # if p1 already in cluster

      clusterm[p2, 2] <- clusterm[p1, 2]

      clust <- which(clusterm[, 2] == clusterm[p1, 2])

      # set distm to infinity to rule out inner cluster distances
      distm[clust, clust] <- Inf

      # update nr clusters
      nr_clusts <- nr_clusts - 1

    }

    # update radius
    all_clust <- unique(clusterm[!is.na(clusterm[,2]),2])

    for(i in all_clust) {

      clust_paths <- which(clusterm[,2] == i)

      if(length(clust_paths) == 2) {

        r <- distm_copy[clust_paths[1], clust_paths[2]]

      } else {

        r <- min(apply(distm_copy[clust_paths, clust_paths], 1, max))

      }

      if(r > rad) rad <- r

    }


  }

  return(clusterm)


}



#### single linkage radius and based ####
radius_k_single <- function(distm, distm_copy, nr_paths, radius, k) {

  # get rid of distances to paths itself - inf

  diag(distm) <- Inf

  # define cluster matrix - col1 = paths, col2 = cluster id

  clusterm <- matrix(0, nr_paths, 2)

  clusterm[, 1] <- 1:nr_paths
  clusterm[, 2] <- NA

  # find clusters based on radius

  #initialize cluster indexing and radius
  nr_clusts <- nr_paths
  clust_counter <- 1
  rad <- 0

  while(rad < radius && nr_clusts > k) {

    min_index <- which(distm == min(distm), arr.ind = TRUE)

    p1 <- min_index[1,1]
    p2 <- min_index[1,2]


    # both not assigned to cluster yet
    if(is.na(clusterm[p1, 2]) && is.na(clusterm[p2, 2])) {

      clusterm[p1, 2] <- clust_counter
      clusterm[p2, 2] <- clust_counter

      distm[p1, p2] <- Inf
      distm[p2, p1] <- Inf

      # update counter nr clusters
      clust_counter <- clust_counter + 1

      # update nr clusters
      nr_clusts <- nr_clusts - 1


    } else if(!is.na(clusterm[p1, 2]) && !is.na(clusterm[p2,2])) {

      # if both already in cluster

      id1 <- clusterm[p1, 2]
      id2 <- clusterm[p2, 2]

      # merge clusters
      clusterm[which(clusterm[,2] == id2), 2] <- id1

      clust <- which(clusterm[,2] == id1)

      # set distm to infinity to rule out inner cluster distances
      distm[clust, clust] <- Inf

      # update nr of clusters
      nr_clusts <- nr_clusts - 1

    } else if(is.na(clusterm[p1, 2])) {

      # if p2 already in cluster

      clusterm[p1, 2] <- clusterm[p2, 2]

      clust <- which(clusterm[,2] == clusterm[p2, 2])

      # set distm to infinity to rule out inner cluster distances
      distm[clust, clust] <- Inf

      # update nr clusters
      nr_clusts <- nr_clusts - 1

    } else {

      # if p1 already in cluster

      clusterm[p2, 2] <- clusterm[p1, 2]

      clust <- which(clusterm[, 2] == clusterm[p1, 2])

      # set distm to infinity to rule out inner cluster distances
      distm[clust, clust] <- Inf

      # update nr clusters
      nr_clusts <- nr_clusts - 1

    }

    # update radius
    all_clust <- unique(clusterm[!is.na(clusterm[,2]),2])

    for(i in all_clust) {

      clust_paths <- which(clusterm[,2] == i)

      if(length(clust_paths) == 2) {

        r <- distm_copy[clust_paths[1], clust_paths[2]]

      } else {

        r <- min(apply(distm_copy[clust_paths, clust_paths], 1, max))

      }

      if(r > rad) rad <- r

    }


  }

  return(clusterm)


}



#### single linkage k based ####
k_single <- function(distm, nr_paths, k) {

  # get rid of distances to paths itself - inf

  diag(distm) <- Inf

  # define cluster matrix - col1 = paths, col2 = cluster id

  clusterm <- matrix(0, nr_paths, 2)

  clusterm[, 1] <- 1:nr_paths
  clusterm[, 2] <- NA

  # find clusters based on k

  #initialize cluster indexing
  nr_clusts <- nr_paths
  clust_counter <- 1

  while(nr_clusts > k) {

    min_index <- which(distm == min(distm), arr.ind = TRUE)

    p1 <- min_index[1,1]
    p2 <- min_index[1,2]


    # both not assigned to cluster yet
    if(is.na(clusterm[p1, 2]) && is.na(clusterm[p2, 2])) {

      clusterm[p1, 2] <- clust_counter
      clusterm[p2, 2] <- clust_counter

      distm[p1, p2] <- Inf
      distm[p2, p1] <- Inf

      # update counter nr clusters
      clust_counter <- clust_counter + 1

      # update nr clusters
      nr_clusts <- nr_clusts - 1


    } else if(!is.na(clusterm[p1, 2]) && !is.na(clusterm[p2,2])) {

      # if both already in cluster

      id1 <- clusterm[p1, 2]
      id2 <- clusterm[p2, 2]

      # merge clusters
      clusterm[which(clusterm[,2] == id2), 2] <- id1

      clust <- which(clusterm[,2] == id1)

      # set distm to infinity to rule out inner cluster distances
      distm[clust, clust] <- Inf

      # update nr of clusters
      nr_clusts <- nr_clusts - 1

    } else if(is.na(clusterm[p1, 2])) {

      # if p2 already in cluster

      clusterm[p1, 2] <- clusterm[p2, 2]

      clust <- which(clusterm[,2] == clusterm[p2, 2])

      # set distm to infinity to rule out inner cluster distances
      distm[clust, clust] <- Inf

      # update nr clusters
      nr_clusts <- nr_clusts - 1

    } else {

      # if p1 already in cluster

      clusterm[p2, 2] <- clusterm[p1, 2]

      clust <- which(clusterm[, 2] == clusterm[p1, 2])

      # set distm to infinity to rule out inner cluster distances
      distm[clust, clust] <- Inf

      # update nr clusters
      nr_clusts <- nr_clusts - 1

    }

  }

  return(clusterm)

}









