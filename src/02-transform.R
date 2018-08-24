source("src/common.R")
source("src/io.R")

library(simIReff)
library(doParallel)
stopImplicitCluster()
registerDoParallel(cores = max(1, detectCores()-1))

for(measure in .MEASURES) {
  for(collection in .COLLECTIONS) {
    cat("\n", collection, measure); flush.console()

    in.path <- file.path("scratch/01-margins/", paste0(collection, "_", measure))
    out.path <- file.path("scratch/02-transform/", paste0(collection, "_", measure))
    dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

    effs.files <- list.files(in.path, full.names = TRUE)
    foreach(ei = seq_along(effs.files), .packages = "simIReff") %dopar% { # for every fitted distribution
      eff.file <- effs.files[ei]
      set.seed(ei)

      eff <- load.object(eff.file) # load and transform
      teff <- try(effTransform(eff, abs.tol = 1e-5), silent = TRUE)
      if(!inherits(teff, "try-error")) {
        save.object(teff, file.path(out.path, basename(eff.file)))
      } else {
        cat(eff.file, " failed.\n")
      }
    }
  }
}

stopImplicitCluster()
