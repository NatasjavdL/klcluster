### curve simplification

#'@title simplify under frechet
#'
#'@description The function simplifies a polygonal curve P by an approximation
#' P' subject to an error epsilon under the Frechet distance.
#'
#'@param path This is a data.frame containing two columns with
#'information about the original trajectory/path. The columns are named x and y
#'and contain the x and y coordinates of the trajectory.
#'
#'@param eps This is a number representing epsilon which will be input for the
#'Frechet distance.
#'
#'@export

simplify_curve <- function(path, eps) {

  size_p <- nrow(path)

  j <- 1
  i <- 1

  new_path <- path[i[j]]

  while(i[j] < size_p) {

    l <- 0

    # add condition in case the index exceeds the number of vertices
    while(2^(l+1) + i[j] <= size_p &
        check_eps_c(as.matrix(path[c(i[j], (i[j] + 2^(l+1))), ]),
          as.matrix(path[i[j]:(i[j] + 2^(l+1)), ]), eps)) {

      l <- l + 1

    }

    low <- 2^l
    high <- 2^(l+1)

    while(low < (high-1)) {

      mid <- (low + high)/2

      if(i[j] + mid > size_p) {

        high <- mid

      } else if(check_eps_c(as.matrix(path[c(i[j], (i[j] + mid)), ]),
        as.matrix(path[i[j]:(i[j] + mid), ]), eps)) {

        low <- mid

      } else {

        high <- mid

      }

    }

    i <- c(i, (i[j] + low))

    j <- j + 1

  }

  new_path <- path[i, ]

  if(sum(new_path[nrow(new_path), ] == path[size_p, ]) != 2) {

    new_path <- rbind(new_path, path[size_p, ])

  }


  return(new_path)

}


