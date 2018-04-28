source("src/common.R")
source("src/io.R")

library(doParallel)
library(rvinecopulib)

stopImplicitCluster()
registerDoParallel(cores = max(1, detectCores()-1))

for(measure in .MEASURES) {
  for(collection in .COLLECTIONS) {
    in.path <- file.path("scratch/06-cop_fit", paste0(collection, "_", measure))
    out.path <- file.path("scratch/07-cop_sim", paste0(collection, "_", measure))
    dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

    dat <- read.csv(file.path("data", paste0(collection, "_", measure, ".csv")))
    dup.i <- indexOfDuplicates(dat)
    if(length(dup.i) > 0)
      dat <- dat[, -dup.i]

    effs <- lapply(colnames(dat), function(f)
      load.object(file.path("scratch/03-margins_select", paste0(collection, "_", measure),
                            "AIC", f)))

    cop <- load.object(file.path(in.path, "AIC_rvine_full"))

    foreach(trial = 1:1000) %dopar% {
      cat("\n", collection, measure, trial); flush.console()

      r <- rvinecop(1000, cop)
      y <- sapply(1:ncol(r), function(i) qeff(r[,i], effs[[i]]))
      colnames(y) <- colnames(dat)
      y <- round(y, 5)
      write.csv(file = file.path(out.path, paste0(trial, ".csv")), row.names = FALSE, y)
    }
  }
}

stopImplicitCluster()
