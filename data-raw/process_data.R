## load data Richard Mann

library(data.table)
library(tidyverse)


pigeon_data <- function(folder, subfolders) {

  for(i in subfolders) {

    new_file_path <- paste0(folder, i)

    files <- list.files(path = new_file_path, pattern = ".txt")

    all_paths <- list()
    all_paths2 <- list()

    for(j in 1:length(files)) {

      dt <- fread(file = paste0(new_file_path, "/", files[j]), fill=TRUE,
        quote='')

      if(names(dt)[1] == "V1") {

        setnames(dt, names(dt), c("Longitude", "Latitude",
          names(dt)[3:length(names(dt))]))

      }

      if(names(dt)[1] != "Longitude") {

        setnames(dt, grep("ongitude", names(dt), value = TRUE), "Longitude")
        setnames(dt, grep("atitude", names(dt), value = TRUE), "Latitude")
      }

      setnames(dt, c("Longitude", "Latitude"), c("long", "lat"))

      all_paths[[j]] <- dt
      all_paths2[[j]] <- dt[, c("long", "lat")]

    }

    assign(i, all_paths2, envir = .GlobalEnv)

  }

}



folder <- paste0(getwd(), "/data-raw/")

subfolders <- c("a55", "brc", "p29", "p94")

pigeon_data(folder, subfolders)

# use data

devtools::use_data(a55)
devtools::use_data(brc)
devtools::use_data(p29)
devtools::use_data(p94)


