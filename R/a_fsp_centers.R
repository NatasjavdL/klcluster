#'@title fsp centers
#'
#'@param paths This is a list of data.frames. Each data.frame represents a path.
#'Each data.frame contains two columns. The columns are named x and y
#'and contain the x and y coordinates of the trajectory
#'
#'@param clust_out the output of one of the clustering algorithms
#'
#'@param dist_measure input either "cont_frechet", "disc_frechet", "dtw or 1,2,3
#'respectively. Specifies the distance measure to be used.
#'
#'@param dtw_settings a vector of 3 variables with settings when using dtw
#'distance c('method_dtw', 'squared_dist' ,' normalize_dist' )
#'method_dtw: integer 1,2,3,4; specified what center type to use
#'squared_dist: bool; specifies wheter to use squared distances
#'normalize_dist: bool; specifies wheteher to normalize distances to output in
#'cost in the same range as frechet distance
#'
#'@return Clustering output containing clusters, centers and result
#'
#'@export

fsp_centers <- function(paths, clust_out, dist_measure = 3,
  dtw_settings = c(1, TRUE, TRUE)) {

  if(is.integer(dist_measure)) {

    tmp <- c("cont_frechet", "disc_frechet", "dtw")
    dist_measure <- tmp[dist_measure]

  }

  if(dist_measure == "dtw") {

    new_out <- fsp_dtw(paths, clust_out, dtw_settings)

  } else if(dist_measure == "disc_frechet") {

    new_out <- fsp_disc_frechet(paths, clust_out)

  } else {

    new_out <- fsp_cont_frechet(paths, clust_out)

  }

  # remove parameters from new_out
  new_out$parameters <- NULL

  # get new distances
  distm <- data.frame("c0" = 1:length(paths))

  for(i in 1:length(new_out$centers)) {

    distm[paste0("c", i)] <- sapply(1:length(paths), function(x)
      get_similarity(dist_measure, new_out$centers[[i]], paths[[x]],
        sqr_dist = dtw_settings[2], norm_dist = dtw_settings[3]))

  }

  distm["c0"] = NULL

  row_mins <- apply(distm, 1, min)

  results <- list()
  results$cost <- max(row_mins)
  results$minsum <- sum(row_mins)
  results$center_ids <- new_out$result$center_ids
  results$distances <- df_to_mat(distm)

  new_out$result <- results

  return(new_out)

}


#@title sample centers
#
#

sample_center <- function(path, size) {

  # initialize matrix
  new_center <- matrix(0, (nrow(path) + size*(nrow(path) - 1)), 2)

  index <- 1

  for(i in 1:(nrow(path) - 1)) {

    # get points defining line segment
    p1 <- path[i, ]
    p2 <- path[i+1, ]

    # sample points
    seg_x <- seq(from = p1[1], to = p2[1], length.out = size + 2)

    seg_y <- seq(from = p1[2], to = p2[2], length.out = size + 2)

    # add to center
    new_index <- index + size + 1
    new_center[index:new_index, 1] <- seg_x
    new_center[index:new_index, 2] <- seg_y

    index <- new_index

  }

  return(new_center)

}

method_dtw <- function(dtw_m, path, center) {

  if(dtw_m == 1) {

    # use original center
    return(center)

  } else if(dtw_m == 2) {

    # use original path
    return(df_to_mat(path))

  } else if(dtw_m == 3) {

    # sample center to equally distribute points to roughly 15% path size
    perc <- nrow(path) * 0.15

    # gives nr points per segment to sample
    size <- ceiling(perc/(nrow(center) - 1))

    return(sample_center(center, size))

  } else if(dtw_m == 4) {

    # simplify center to 10% path size
    perc <- ceiling(nrow(path) * 0.10)

    return(df_to_mat(bin_simplify(path, perc, 1)))

  }

}


fsp_dtw <- function(paths, clust_out, dist_settings) {

  # get input from clust_out
  centers <- clust_out$centers
  clusters <- clust_out$clusters
  ids <- clust_out$result$center_ids
  method <- dist_settings[1]

  new_centers <- list()

  for(i in 1:length(clusters)) {

    # subset paths and centers
    clust <- clusters[[i]]
    paths_used <- paths[clust]

    # skip iteration when size cluster = 1 - only one center possible
    if(length(clust) == 1) {

      new_centers[[i]] <- centers[[i]]

      next

    }

    # get center to use for optimizing quality - depends on method
    c <- method_dtw(method, paths[[ids[i]]], df_to_mat(centers[[i]]))

    # initialize points matrix
    pm <- matrix(0, 0, 3)

    # gather points from distance matrix
    for(j in 1:length(clust)) {

      # select path
      p <- df_to_mat(paths_used[[j]])

      # get distance matrix
      dm <- dtw(c, p, dist_settings[2])

      # backtrack points from dm
      points <- backtrack(c, p, dm)

      # concatenate points to points matrix
      pm <- rbind(pm, points)
    }

    # initialize new center matrix
    c_new <- matrix(0, nrow(c), 2)

    # for each coordinate of center, get mean of points aligned with coordinate
    for(k in 1:nrow(c)) {

      subset_pm <- pm[(pm[,1] == k), 2:3]

      # if only one mapped - assign directly
      if(sum(pm[,1] == k) == 1) {

        c_new[k, ] <- subset_pm

        # else calculate coordinate wise mean
      } else {

        c_new[k, ] <- apply(subset_pm, 2, sum)/nrow(subset_pm)

      }

    }

    # change new center to data frame
    c_new <- mat_to_df(c_new, path = TRUE)


    # simplify center to 10% of center if original is used and method = 2
    if(method == 2) {

      size <- floor(nrow(c_new) * 0.05)
      c_new <- bin_simplify(c_new, size, 1)

    }

    new_centers[[i]] <- c_new

  }

  # update clust_out
  clust_out$centers <- new_centers

  return(clust_out)

}


fsp_disc_frechet <- function(paths, clust_out) {

  # get input from clust_out
  centers <- clust_out$centers
  clusters <- clust_out$clusters

  new_centers <- list()

  for(i in 1:length(clusters)) {

    # subset paths and centers
    clust <- clusters[[i]]
    paths_used <- paths[clust]

    # if only one path in cluster, skip and keep original center
    if(length(clust) == 1) {

      new_centers[[i]] <- centers[[i]]

      # skip iteration
      next

    }

    c <- df_to_mat(centers[[i]])

    # initialize points matrix
    pm <- matrix(0, 0, 3)

    # gather points from distance matrix
    for(j in 1:length(clust)) {

      # select path
      p <- df_to_mat(paths_used[[j]])

      # get distance matrix
      dm <- dist_matrix(c, p)

      # backtrack points from dm
      points <- backtrack(c, p, dm)

      # concatenate points to points matrix
      pm <- rbind(pm, points)
    }

    c_new <- matrix(0, nrow(c), 2)

    for(k in 1:nrow(c)) {

      subset_pm <- pm[(pm[,1] == k), 2:3]

      if(sum((pm[,1] == k)) == 1) {

        c_new[k, ] <- subset_pm

      } else {

        c_new[k, ] <- sed(subset_pm)

      }

    }

    c_new <- mat_to_df(c_new, path = TRUE)
    new_centers[[i]] <- c_new

  }

  clust_out$centers <- new_centers

  return(clust_out)

}

fsp_cont_frechet <- function(paths, clust_out) {

  centers <- clust_out$centers
  clusters <- clust_out$clusters

  new_centers <- list()

  for(i in 1:length(clusters)) {

    # subset paths and centers
    clust <- clusters[[i]]

    paths_used <- paths[clust]


    if(length(clust) == 1) {

      new_centers[[i]] <- centers[[i]]

      # skip iteration
      next

    }

    c <- df_to_mat(centers[[i]])

    # initialize points matrix
    pm <- matrix(0, 0, 3)

    # gather points from distance matrix
    for(j in 1:length(clust)) {

      # select path
      p <- df_to_mat(paths_used[[j]])

      # get eps value (plus 1.1 to account for precision)
      eps <- frechet_distance(c, p, 0.1) + 0.15

      # get coordinates
      points <- segment_wise_coords(c, p, eps)

      # concatenate points to points matrix
      pm <- rbind(pm, points)
    }

    c_new <- matrix(0, nrow(c), 2)

    for(k in 1:nrow(c)) {

      subset_pm <- pm[(pm[,1] == k), 2:3]

      c_new[k, ] <- sed(subset_pm)

    }

    c_new <- mat_to_df(c_new, path = TRUE)
    new_centers[[i]] <- c_new

  }

  clust_out$centers <- new_centers

  return(clust_out)


}

