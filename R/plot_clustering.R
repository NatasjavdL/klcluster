## functions for plotting clusterings

#'@title plot all clusters without dots
#'
#'@param all_paths is the list of data.frames with the trajectories
#'
#'@param out_dt contains the clustering output
#'
#'@param name_plot is a string containing the title of the plot
#'
#'@export

plot_grid_no_points <- function(all_paths, out_dt, name_plot = "") {

  paths <- as.data.frame(data.table::rbindlist(all_paths, idcol = TRUE))
  nr_clusts <- length(out_dt$clusters)
  paths$clust <- 0

  for(j in 1:nr_clusts) {

    cluster <- out_dt$clusters[[j]]

    paths$clust[(paths$`.id` %in% cluster)] <- j

  }

  counter <- 1
  fp <- data.frame(matrix(0, 0, 3))
  names(fp) <- c("long", "lat", "id")
  clusts <- out_dt$clusters

  for(j in 1:nr_clusts) {

    paths_in_clust <- clusts[[j]]

    for(i in paths_in_clust) {

      fp[counter, 1:2] <- all_paths[[i]][1,]
      fp[counter, 3] <- j
      counter <- counter + 1

    }

  }

  col_pallete <- colorspace::rainbow_hcl(nr_clusts)

  plot_clust <- ggplot()

  grid_plot <- list()

  for(i in 1:nr_clusts) {

    paths_used <- paths[(paths$clust == i), ]

    path_plot <- geom_path(data = paths_used, alpha = 0.8,
      aes(x = long, y = lat, group = .id, col = clust), color = col_pallete[i])
    center_plot <- geom_path(data = out_dt$centers[[i]],
        aes(x = long, y = lat), col = "black")
    center_points <- geom_point(data = out_dt$centers[[i]],
        aes(x = long, y = lat), col = "black", shape = ".")
    fp_plot <- geom_point(data = fp[fp$id == i,], aes(x = long, y = lat),
        shape = 23, fill="red", color="black", size=2)

    plot_clust <- plot_clust + path_plot + center_plot + center_points +
    fp_plot

    grid_plot[[i]] <- ggplot() + path_plot + center_plot + center_points +
      fp_plot +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

  }

  plot_clust <- plot_clust +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title=element_blank(),
      axis.text=element_blank(),
      axis.ticks=element_blank()) +
    ggtitle(name_plot)

  grid1 <- cowplot::plot_grid(plotlist = grid_plot, align = 'h', nrow = 1)
  final_grid <- cowplot::plot_grid(grid1, plot_clust, align = 'v', nrow = 2,
    rel_heights = c(1, 3))

  return(final_grid)

}

#'@title plot all clusters with points and grid
#'
#'@param all_paths is the list of data.frames with the trajectories
#'
#'@param out_dt contains the clustering output
#'
#'@param name_plot is a string containing the title of the plot
#'
#'@export

plot_all_clust_grid <- function(all_paths, out_dt, name_plot = "") {

  paths <- as.data.frame(data.table::rbindlist(all_paths, idcol = TRUE))
  nr_clusts <- length(out_dt$clusters)
  paths$clust <- 0

  for(j in 1:nr_clusts) {

    cluster <- out_dt$clusters[[j]]

    paths$clust[(paths$`.id` %in% cluster)] <- j

  }

  counter <- 1
  fp <- data.frame(matrix(0, 0, 3))
  names(fp) <- c("long", "lat", "id")
  clusts <- out_dt$clusters

  for(j in 1:nr_clusts) {

    paths_in_clust <- clusts[[j]]

    for(i in paths_in_clust) {

      fp[counter, 1:2] <- all_paths[[i]][1,]
      fp[counter, 3] <- j
      counter <- counter + 1

    }

  }

  col_pallete <- colorspace::rainbow_hcl(nr_clusts)

  plot_clust <- ggplot()

  grid_plot <- list()

  for(i in 1:nr_clusts) {

    paths_used <- paths[(paths$clust == i), ]

    path_plot <- geom_path(data = paths_used, alpha = 0.8,
      aes(x = long, y = lat, group = .id, col = clust), color = col_pallete[i])
    center_plot <- geom_path(data = out_dt$centers[[i]],
      aes(x = long, y = lat), col = "black")
    center_points <- geom_point(data = out_dt$centers[[i]],
      aes(x = long, y = lat), shape = 21, fill = "black", color = col_pallete[i])
    fp_plot <- geom_point(data = fp[fp$id == i,], aes(x = long, y = lat),
      shape = 23, fill="red", color="black", size=2)

    plot_clust <- plot_clust + path_plot + center_plot + center_points +
      fp_plot

    grid_plot[[i]] <- ggplot() + path_plot + center_plot + center_points +
      fp_plot +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

  }

  plot_clust <- plot_clust +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title=element_blank(),
      axis.text=element_blank(),
      axis.ticks=element_blank()) +
    ggtitle(name_plot)

  grid1 <- cowplot::plot_grid(plotlist = grid_plot, align = 'h', nrow = 1)
  final_grid <- cowplot::plot_grid(grid1, plot_clust, align = 'v', nrow = 2,
    rel_heights = c(1, 3))

  return(final_grid)

}


#'@title plot all clusters without grid
#'
#'@param all_paths is the list of data.frames with the trajectories
#'
#'@param out_dt contains the clustering output
#'
#'@param name_plot is a string containing the title of the plot
#'
#'@export

plot_all_clust <- function(all_paths, out_dt, name_plot = "") {

  paths <- as.data.frame(data.table::rbindlist(all_paths, idcol = TRUE))
  nr_clusts <- length(out_dt$clusters)
  paths$clust <- 0

  for(j in 1:nr_clusts) {

    cluster <- out_dt$clusters[[j]]

    paths$clust[(paths$`.id` %in% cluster)] <- j

  }

  counter <- 1
  fp <- data.frame(matrix(0, 0, 3))
  names(fp) <- c("long", "lat", "id")
  clusts <- out_dt$clusters

  for(j in 1:nr_clusts) {

    paths_in_clust <- clusts[[j]]

    for(i in paths_in_clust) {

      fp[counter, 1:2] <- all_paths[[i]][1,]
      fp[counter, 3] <- j
      counter <- counter + 1

    }

  }

  col_pallete <- colorspace::rainbow_hcl(nr_clusts)

  plot_clust <- ggplot()

  for(i in 1:nr_clusts) {

    paths_used <- paths[(paths$clust == i), ]

    path_plot <- geom_path(data = paths_used,
      aes(x = long, y = lat, group = .id, col = clust), color = col_pallete[i])
    center_plot <- geom_path(data = out_dt$centers[[i]],
      aes(x = long, y = lat), col = "black")
    center_points <- geom_point(data = out_dt$centers[[i]],
      aes(x = long, y = lat), col = "black")
    fp_plot <- geom_point(data = fp[fp$id == i,], aes(x = long, y = lat),
      shape = 23, fill="red", color="black", size=3)

    plot_clust <- plot_clust + path_plot + center_plot + center_points + fp_plot
  }

  plot_clust <- plot_clust +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(name_plot)

  return(plot_clust)

}



#'@title plot all clusters on a mao using ggmap::ggmap and stamenmaps
#'
#'@param all_paths is the list of data.frames with the trajectories
#'
#'@param out_dt contains the clustering output
#'
#'@param name_plot is a string containing the title of the plot
#'
#'@param map_type contains the map type you want to use, see
#'[stamenmaps](maps.stamen.com)
#'
#'@param zoom_ggmap defines the quality of the map
#'
#'@export

plot_all_clust_map <- function(all_paths, out_dt, name_plot = "",
  map_type, zoom_ggmap = 13) {

  # need to get a bounding box - use any path and extend if necessary
  tmp <- do.call(rbind, all_paths)
  # names(tmp) <- c("x", "y") # change names long/lat to x/y for extent
  # m <- get_map(bbox(extent(tmp)*factor_bbox), zoom=zoom_ggmap)

  height <- max(tmp$lat) - min(tmp$lat)
  width <- max(tmp$long) - min(tmp$long)
  map_borders <- c(bottom  = min(tmp$lat)  - 0.1 * height,
    top     = max(tmp$lat)  + 0.1 * height,
    left    = min(tmp$long) - 0.1 * width,
    right   = max(tmp$long) + 0.1 * width)
  map <- ggmap::get_stamenmap(map_borders, maptype = map_type, zoom = zoom_ggmap)

  # transform dataset
  paths <- as.data.frame(data.table::rbindlist(all_paths, idcol = TRUE))
  nr_clusts <- length(out_dt$clusters)
  paths$clust <- 0

  for(j in 1:nr_clusts) {

    cluster <- out_dt$clusters[[j]]

    paths$clust[(paths$`.id` %in% cluster)] <- j

  }

  counter <- 1
  fp <- data.frame(matrix(0, 0, 3))
  names(fp) <- c("long", "lat", "id")
  clusts <- out_dt$clusters

  for(j in 1:nr_clusts) {

    paths_in_clust <- clusts[[j]]

    for(i in paths_in_clust) {

      fp[counter, 1:2] <- all_paths[[i]][1,]
      fp[counter, 3] <- j
      counter <- counter + 1

    }

  }

  # plot clustering
  col_pallete <- colorspace::rainbow_hcl(nr_clusts)

  # plot_clust <- ggmap::ggmap(m)
  plot_clust <- ggmap::ggmap(map)

  for(i in 1:nr_clusts) {

    paths_used <- paths[(paths$clust == i), ]

    plot_clust <- plot_clust + geom_path(data = paths_used,
      aes(x = long, y = lat, group = .id, col = clust), color = col_pallete[i]) +
      geom_path(data = out_dt$centers[[i]],
        aes(x = long, y = lat), col = "black") +
      geom_point(data = out_dt$centers[[i]],
        aes(x = long, y = lat), shape = 21, fill = "black", color = col_pallete[i]) +
      geom_point(data = fp[fp$id == i,], aes(x = long, y = lat),
        shape = 23, fill="red", color="black", size=3)
  }

  plot_clust <- plot_clust +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(name_plot)

  return(plot_clust)

}



#'@title plot all centers on maps
#'
#'@param all_paths is the list of data.frames with the trajectories
#'
#'@param out_dt contains the clustering output
#'
#'@param name_plot is a string containing the title of the plot
#'
#'@param map_type contains the map type you want to use, see
#'[stamenmaps](maps.stamen.com)
#'
#'@param zoom_ggmap defines the quality of the map
#'
#'@export

plot_centers_map <- function(all_paths, out_dt, name_plot = "",
  map_type = "terrain", zoom_ggmap = 13) {

  # need to get a bounding box - use any path and extend if necessary
  tmp <- do.call(rbind, all_paths)
  # names(tmp) <- c("x", "y") # change names long/lat to x/y for extent
  # m <- get_map(bbox(extent(tmp)*factor_bbox), zoom=zoom_ggmap)

  height <- max(tmp$lat) - min(tmp$lat)
  width <- max(tmp$long) - min(tmp$long)
  map_borders <- c(bottom  = min(tmp$lat)  - 0.1 * height,
    top     = max(tmp$lat)  + 0.1 * height,
    left    = min(tmp$long) - 0.1 * width,
    right   = max(tmp$long) + 0.1 * width)
  map <- ggmap::get_stamenmap(map_borders, maptype = map_type, zoom = zoom_ggmap,
    messaging = FALSE)

  # transform dataset
  paths <- as.data.frame(data.table::rbindlist(all_paths, idcol = TRUE))
  nr_clusts <- length(out_dt$clusters)
  paths$clust <- 0

  for(j in 1:nr_clusts) {

    cluster <- out_dt$clusters[[j]]

    paths$clust[(paths$`.id` %in% cluster)] <- j

  }

  counter <- 1
  fp <- data.frame(matrix(0, 0, 3))
  names(fp) <- c("long", "lat", "id")
  clusts <- out_dt$clusters

  for(j in 1:nr_clusts) {

    paths_in_clust <- clusts[[j]]

    for(i in paths_in_clust) {

      fp[counter, 1:2] <- all_paths[[i]][1,]
      fp[counter, 3] <- j
      counter <- counter + 1

    }

  }

  # plot clustering
  col_pallete <- colorspace::rainbow_hcl(nr_clusts)

  # plot_clust <- ggmap::ggmap(m)
  plot_clust <- ggmap::ggmap(map)

  for(i in 1:nr_clusts) {

    paths_used <- paths[(paths$clust == i), ]

    plot_clust <- plot_clust +
      geom_path(data = out_dt$centers[[i]],
        aes(x = long, y = lat), col = "black") +
      geom_point(data = out_dt$centers[[i]], aes(x = long, y = lat), shape = 21,
        color = "black", fill = col_pallete[i], size = 3) +
      geom_point(data = fp[fp$id == i,], aes(x = long, y = lat),
        shape = 23, fill="red", color="black", size=3)
  }

  plot_clust <- plot_clust +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ggtitle(name_plot)

  return(plot_clust)

}

