source("src/common.R")
source("src/io.R")

library(parallel)
library(VineCopula)
library(rvinecopulib)

for(measure in .MEASURES) {
  for(collection in .COLLECTIONS) {
    dat <- read.csv(file.path("data", paste0(collection, "_", measure, ".csv"))) # original data
    dup.i <- indexOfDuplicates(dat) # remove identical systems
    if(length(dup.i) > 0)
      dat <- dat[, -dup.i]

    set.seed(ncol(dat))
    U <- pobs(dat, ties.method = "random") # pseudo-obsevations from original data

    out.path <- file.path("scratch/06-cop_fit/", paste0(collection, "_", measure))
    dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

    for(selcrit in c("AIC", "BIC", "logLik")) { # for each selection criteria
      cat(selcrit, collection, measure, "\n"); flush.console()

      cop <- vinecop(U, family_set = "gauss", selcrit = tolower(selcrit), cores = 1, presel = FALSE)
      save.object(cop, file.path(out.path, paste0(selcrit, "_gaussian_full"))) # gaussian

      cop <- vinecop(U, family_set = c("onepar", "twopar"), selcrit = tolower(selcrit),
                     cores = max(1, detectCores()-1), presel = FALSE) # r-vine
      save.object(cop, file.path(out.path, paste0(selcrit, "_rvine_full")))
    }
  }
}
