source("src/common.R")
source("src/io.R")

library(rvinecopulib)
library(doParallel)
stopImplicitCluster()
registerDoParallel(cores = max(1, detectCores()-1))

for(measure in .MEASURES) {
  for(collection in .COLLECTIONS) {
    in.path <- file.path("scratch/06-cop_fit", paste0(collection, "_", measure))
    out.path <- file.path("output/07-cop_sim")
    dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

    dat <- read.csv(file.path("data", paste0(collection, "_", measure, ".csv"))) # original data
    dup.i <- indexOfDuplicates(dat) # remove identical systems
    if(length(dup.i) > 0)
      dat <- dat[, -dup.i]

    effs <- lapply(colnames(dat), function(f) # load all margins
      load.object(file.path("scratch/03-margins_select", paste0(collection, "_", measure),
                            "AIC", f)))

    cop <- load.object(file.path(in.path, "AIC_rvine_full")) # load copula

    r <- foreach(trial = 1:1000, .combine = rbind, .inorder = TRUE) %dopar% {
      set.seed(trial)
      cat("\n", collection, measure, trial); flush.console()

      r <- rvinecop(1000, cop) # simulate pseudo-observations
      y <- sapply(1:ncol(r), function(i) qeff(r[,i], effs[[i]])) # pass through quantile functions
      colnames(y) <- colnames(dat)
      apply(y, 2, function(yy) c(mean(yy), var(yy))) # compute per-system mean and variance
    }
    
    r <- round(r, 8)
    write.csv(file = file.path(out.path, paste0(collection, "_", measure, "_mean.csv")),
              row.names = FALSE, r[seq(1, nrow(r), 2),])
    write.csv(file = file.path(out.path, paste0(collection, "_", measure, "_var.csv")),
              row.names = FALSE, r[seq(2, nrow(r), 2),])
  }
}

stopImplicitCluster()
