#'@title Transform movebank data to input format
#'
#'@description this function transforms movebank data format to the right input
#'format for the clustering algorithms
#'
#'@param move_dt move data structure, e.g. "Move", "MoveStack", "MoveBurst"
#'
#'@export
#'


transform_move <- function(move_dt) {

  dt_type <- class(move_dt)[1]

  if(!(dt_type %in% c("Move", "MoveStack", "MoveBurst"))) {

    message("data is not in movebank format. More information can be found on
      the following website:
      https://www.rdocumentation.org/packages/move/versions/3.0.1/topics/move")

    return(move_dt)

  }

  if(dt_type == "Move") {

    message("Your data only consists of one path. For clustering please stack
      the paths by using function moveStack(), see
      https://www.rdocumentation.org/packages/move/versions/3.0.1/topics/move")

    return(list(as.data.frame(move_dt)[, c("location_lat", "location_long")]))

  } else {

    # unstack paths
    tmp <- move::split(x = move_dt)

    all_paths <- list()

    for(i in 1:length(tmp)) {

      tmp2 <- as.data.frame(tmp[[i]])[, c("location_lat", "location_long")]

      tmp2 <- setNames(tmp2, c("lat", "long"))

      all_paths[[i]] <- tmp2

    }

    return(all_paths)

  }

}
