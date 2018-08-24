source("src/common.R")
source("src/io.R")

library(simIReff)

for(measure in .MEASURES) {
  for(collection in .COLLECTIONS) {
    cat("\n", collection, measure); flush.console()

    in.path <- file.path("scratch/01-margins/", paste0(collection, "_", measure))
    out1.path <- file.path("scratch/03-margins_select/", paste0(collection, "_", measure))
    out2.path <- file.path("output/03-margins_select/")
    dir.create(out1.path, recursive = TRUE, showWarnings = FALSE)
    dir.create(out2.path, recursive = TRUE, showWarnings = FALSE)

    dat <- read.csv(paste0("data/", collection, "_", measure, ".csv")) # original data

    fits <- matrix("", ncol = 3, nrow = ncol(dat), dimnames = list(NULL, c("logLik", "AIC", "BIC")))
    for(i in 1:nrow(fits)) { # for every system
      sysname <- colnames(dat)[i]
      set.seed(i)

      # load all distributions fitted for the system and select
      eff.files <- list.files(in.path, pattern = paste0("^", sysname, "_.+$"), full.names = TRUE)
      effs <- lapply(eff.files, load.object)
      effs <- sample(effs) # shuffle to minimize biases due to order
      for(method in colnames(fits)) { # for every selection criteria
        eff <- effSelect(effs, method = method)
        fits[i, method] <- eff$model$type

        save.object(eff, file.path(out1.path, method, sysname)) # save selected
      }
    }

    write.csv(file = file.path(out2.path, paste0(collection, "_", measure, ".csv")),
              row.names = FALSE, fits)
  }
}
