#'@title compute 1-NN clusters
#'
#'@param dt list of data frames, each data frame representing a path
#'
#'@param l number of datapoints to return for each center
#'
#'@param prec_l precision on l, i.e. l - prec_l <= |center| <= l + prec_l
#'
#'@param dist_measure either input 1, 2 or 3 correspondig to
#'"continuous frechet", "discrete frechet" or "dtw"
#'
#'@param prec_eps if dist_measure = 1 (continuous frechet) then this is the
#'precision of the epsilon approximation
#'default is NA
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


kNN_trajectories <- function(dt, l, prec_l = 1, dist_measure = 2, prec_eps = NA,
  square_dist = FALSE, norm_dist = TRUE) {

  nr_paths <- length(dt)

  distm <- matrix(NA_integer_, nr_paths, nr_paths)

  for(i in 1:(nr_paths - 1)) {

    j <- i + 1

    tmp <- sapply(j:nr_paths, function(x) get_similarity(dist_measure, dt[[i]],
      dt[[x]], square_dist, norm_dist))

    distm[j:nr_paths, i] <- tmp
    distm[i, j:nr_paths] <- tmp

  }

  diag(distm) <- Inf

  tmp <- apply(distm, 1, min)
  NN_vec <- sapply(1:nr_paths, function(x) which(distm[,x] == tmp[x]))
  id_vec <- 1:nr_paths
  clust_vec <- rep(NA_integer_, nr_paths)

  nr_clust <- 1

  for(i in 1:nr_paths) {

    clust_p <- clust_vec[i]
    clust_NN <- clust_vec[NN_vec[i]]

    # if p_i not assigned to cluster
    if(is.na(clust_p)) {

      if(is.na(clust_NN)) {

        # if p_i & NN_i not assigned to cluster - add p_i and NN_i to new cluster
        clust_vec[i] <- nr_clust
        clust_vec[NN_vec[i]] <- nr_clust

        nr_clust <- nr_clust + 1

      } else {

        # if NN_i contained in c_k, add p_i to c_k
        clust_vec[i] <- clust_NN

      }

    } else {

      if(is.na(clust_NN)) {

        # if p_i assigned to cluster c_k, NN_i unassigned, asssign NN_i to c_k
        clust_vec[NN_vec[i]] <- clust_p

      } else {

        # p_i assigned to c_k, NN_i to c_j, merge clusters

        if(clust_NN != clust_p) {

          tmp <- which(clust_vec == clust_NN)
          clust_vec[tmp] <- clust_p

        }

      }

    }

  }


  # adjust cluster numbering if not sequential
  clust_ids <- unique(clust_vec)
  max_clust <- max(clust_ids)
  if(max_clust > length(clust_ids)) {

    tmp <- rep(NA_integer_, max_clust)
    tmp[clust_ids] <- 1:length(clust_ids)
    clust_vec_new <- mapply(function(x) x = tmp[x], x = clust_vec)
    clust_vec <- clust_vec_new

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

      #subset distance matrix
      dm <- distm[cluster, cluster]

      diag(dm) <- 0

      # use minimum of sum of distaces to determine center
      center_id <- which(colSums(dm) == min(colSums(dm)))[1]

      path_id <- cluster[center_id]

    }

    path_ids_all[i] <- path_id
    center <- dt[[path_id]]
    center <- bin_simplify(center, l, prec_l)
    center_paths[[i]] <- center

    center_dist[, i] <- sapply(1:nr_paths, function(x)
      get_similarity(dist_measure, center, dt[[x]], prec_eps, square_dist,
        norm_dist))

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
  out$clusters <- clust_list
  out$result <- results
  out$parameters <- params


  # return center_paths

  return(out)

}
