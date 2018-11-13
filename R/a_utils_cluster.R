sample_prob <- function(x, seed = TRUE) {

  if(seed) {

    set.seed(1234)

  }

  sample(1:length(x), 1, prob = x/sum(x), replace = T)

}


min_max <- function(x) {

  max_id <- which(x == max(x))

  if(length(max_id) > 1) {

    max_id <- sample(max_id, 1)

  }

  return(max_id)

}



### cluster discrete ###

#@title cluster discrete
#
#
get_cluster_d <- function(centers, paths) {

  nr_paths <- length(paths)
  nr_centers <- length(centers)

  clusters <- list()

  if(nr_centers == 1) {

    clusters[[1]] = 1:length(paths)

    return(clusters)

  }

  dist_m <- matrix(0, nr_paths, nr_centers)

  for(i in 1:nr_centers) {

    dist_m[, i] <- sapply(1:nr_paths, function(x)
      discrete_frechet(centers[[i]], paths[[x]]))

  }

  mincol <- apply(dist_m, 1, function(x) min(x))

  for(i in 1:nr_centers) {

    center_path <- which(dist_m[, i] == min(dist_m[, i]))[1]
    cluster_paths <- which(dist_m[, i] == mincol)

    clusters[[i]] <- cluster_paths

  }

  return(clusters)

}

### cluster continuous ###
#@title cluster continuous
#
#

get_cluster_c <- function(centers, paths, prec = 1) {

  nr_paths <- length(paths)
  nr_centers <- length(centers)

  clusters <- list()

  if(nr_centers == 1) {

    clusters[[1]] = 1:length(paths)

    return(clusters)

  }

  dist_m <- matrix(0, nr_paths, nr_centers)

  for(i in 1:nr_centers) {

    dist_m[, i] <- sapply(1:nr_paths, function(x)
      get_frechet_distance(centers[[i]], paths[[x]], prec))

  }

  mincol <- apply(dist_m, 1, function(x) min(x))

  for(i in 1:nr_centers) {

    center_path <- which(dist_m[, i] == min(dist_m[, i]))[1]
    cluster_paths <- which(dist_m[, i] == mincol)

    clusters[[i]] <- cluster_paths

  }

  return(clusters)

}


### get clusters from distance matrix
#@title cluster from distance matrix



get_clusters <- function(dist_m) {

  row_mins <- apply(dist_m, 1, min)

  clusters <- list()

  for(i in 1:ncol(dist_m)) {

    clusters[[i]] <- which(row_mins == dist_m[, i])

  }

  names(clusters) <- paste0("cluster_", 1:ncol(dist_m))

  return(clusters)
}

