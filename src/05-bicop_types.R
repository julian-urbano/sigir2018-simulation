source("src/common.R")
source("src/io.R")

library(simIReff)
library(VineCopula)
library(doParallel)
stopImplicitCluster()
registerDoParallel(cores = max(1, detectCores()-1))

out.path <- "output/05-bicop_types"
dir.create(out.path, showWarnings = FALSE, recursive = TRUE)

for(measure in .MEASURES) {
  foreach(ci = seq_along(.COLLECTIONS), .packages = c("simIReff", "VineCopula")) %dopar% {
    collection <- .COLLECTIONS[ci]
    cat("\n", collection, measure); flush.console()

    in.path <- file.path("scratch/03-margins_select/", paste0(collection, "_", measure), "AIC")

    dat <- read.csv(file.path("data", paste0(collection, "_", measure, ".csv"))) # original data
    dup.i <- indexOfDuplicates(dat) # remove identical systems
    if(length(dup.i) > 0)
      dat <- dat[, -dup.i]

    cop.types <- data.frame(si = character(0), sj = character(0),
                            logLik = numeric(0), AIC = numeric(0), BIC = numeric(0),
                            stringsAsFactors = FALSE)
    for(i in 1:(ncol(dat)-1)) { # for every sysem si
      set.seed(i) # for reproducibility
      si <- colnames(dat)[i]
      xi <- dat[,i]
      ui <- pobs(xi, ties.method = "random")

      for(j in (i+1):ncol(dat)) { # for every other system sj
        sj <- colnames(dat)[j]
        xj <- dat[,j]
        uj <- pobs(xj, ties.method = "random")

        # fit and select bicopulas by different criteria
        b.logLik <- BiCopSelect(ui, uj, selectioncrit = "logLik")
        b.AIC <- BiCopSelect(ui, uj, selectioncrit = "AIC")
        b.BIC <- BiCopSelect(ui, uj, selectioncrit = "BIC")
        cop.types <- rbind(cop.types, list(si = si, sj = sj,
                                           logLik = b.logLik$family, AIC = b.AIC$family,
                                           BIC = b.BIC$family),
                           stringsAsFactors = FALSE)
      }
    }
    write.csv(file = file.path(out.path, paste0(collection, "_", measure, ".csv")),
              row.names = FALSE, cop.types)
  }
}

stopImplicitCluster()
