source("src/common.R")
source("src/io.R")

library(doParallel)
stopImplicitCluster()
registerDoParallel(cores = max(1, detectCores()-1))

for(measure in .MEASURES) {
  foreach(ci = seq_along(.COLLECTIONS)) %dopar% {
    collection <- .COLLECTIONS[ci]
    cat("\n", collection, measure); flush.console()
    set.seed(ci)

    out.path <- file.path("scratch/01-margins/", paste0(collection, "_", measure))
    dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

    dat <- read.csv(paste0("data/", collection, "_", measure, ".csv")) # original data

    for(i in 1:ncol(dat)) { # for every system
      x <- dat[,i]

      if(measure %in% c("ap", "ndcg20", "err20")) {
        effs <- effContFit(x)
      }else{
        s <- support(measure)
        effs <- effDiscFit(x, s)
      }

      for(eff in effs) { # save all fitted distributions
        save.object(eff, file.path(out.path, paste0(colnames(dat)[i], "_", eff$model$type)))
      }
    }
  }
}

stopImplicitCluster()
