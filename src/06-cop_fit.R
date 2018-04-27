source("src/common.R")
source("src/io.R")

library(parallel)
library(VineCopula)
library(rvinecopulib)

set.seed(1)
for(measure in .MEASURES) {
  for(collection in .COLLECTIONS) {
    dat <- read.csv(file.path("data", paste0(collection, "_", measure, ".csv")))
    dup.i <- indexOfDuplicates(dat)
    if(length(dup.i) > 0)
      dat <- dat[, -dup.i]

    U <- pobs(dat, ties.method = "random")

    out.path <- file.path("scratch/06-cop_fit/", paste0(collection, "_", measure))
    dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

    for(selcrit in c("AIC", "BIC", "logLik")) {
      cat(selcrit, collection, measure, "\n"); flush.console()

      cop <- vinecop(U, family_set = "gauss", selcrit = tolower(selcrit), cores = 1, presel = FALSE)
      #RVineStructureSelect(data = U, cores = 1, selectioncrit = selcrit, familyset = 1)
      save.object(cop, file.path(out.path, paste0(selcrit, "_gaussian_full")))

      cop <- vinecop(U, family_set = c("onepar", "twopar"), selcrit = tolower(selcrit),
                     cores = max(1, etectCores()-1), presel = FALSE)
      #cop <- RVineStructureSelect(data = U, cores = detectCores()-1, selectioncrit = selcrit)
      save.object(cop, file.path(out.path, paste0(selcrit, "_rvine_full")))
      # cop <- RVineStructureSelect(data = U, cores = detectCores()-1,
      #selectioncrit = selcrit, trunclevel = 5)
      # save.object(cop, file.path(out.path, paste0(selcrit, "_rvine_5")))
    }
  }
}
