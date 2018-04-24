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

    dat <- read.csv(paste0("data/", collection, "_", measure, ".csv"))

    for(i in 1:ncol(dat)) {
      x <- dat[,i]

      if(measure %in% c("ap", "ndcg20", "err20")) {
        effs <- effContFit(x)
      }else{
        support <- switch (measure,
                           "p10" = seq(0, 1, .1),
                           "p20" = seq(0, 1, .05),
                           "rr" = c(0, 1/1000:1)
        );
        effs <- effDiscFit(x, support)
      }

      for(eff in effs) {
        save.object(eff, file.path(out.path, paste0(colnames(dat)[i], "_", eff$model$type)))
      }
    }
  }
}

stopImplicitCluster()
