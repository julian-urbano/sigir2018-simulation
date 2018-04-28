source("src/common.R")
source("src/io.R")

library(doParallel)

out.path <- file.path("output/08-cop_sim_deltas")
dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

stopImplicitCluster()
registerDoParallel(cores = max(1, detectCores()-1))

for(collection in .COLLECTIONS) {
  for(measure in .MEASURES) {
    cat("\n", collection, measure); flush.console()
    # measure <- "ap"
    # collection <- "web2014"

    in.path.eff <- file.path("scratch/03-margins_select/", paste0(collection, "_", measure), "AIC")
    in.path.sim <- file.path("scratch/07-cop_sim/", paste0(collection, "_", measure))

    sysnames <- colnames(read.csv(file.path(in.path.sim, "1.csv")))
    effs <- lapply(sysnames, function(sysname)
      load.object(file.path(in.path.eff, sysname)))

    E <- sapply(effs, function(e) e$mean)
    Var <- sapply(effs, function(e) e$var)

    r <- foreach(i = 1:1000, .combine = "rbind") %dopar% {
      dat <- read.csv(file.path(in.path.sim, paste0(i, ".csv")))
      unname(t(apply(dat, 2, function(d) c(mean(d), var(d)))))
    }
    r <- data.frame(sys = seq_along(effs), mean = r[,1]-E, var = r[,2]-Var)

    write.csv(file = file.path(out.path, paste0(collection, "_", measure, ".csv")),
              row.names = FALSE, r)
  }
}
q(save = "no")

n.trials <- 1000
n.sample <- 1000
deltas <- replicate(n.trials, {
  shape1 <- runif(1, 1, 10)
  shape2 <- runif(1, 1, 10)
  E <- shape1 / (shape1 + shape2)
  Var <- shape1 * shape2 / (shape1 + shape2)^2 / (shape1 + shape2 + 1)
  x <- rbeta(n.sample, shape1, shape2)
  c(mean(x) - E, var(x) - Var)
})

deltas <- replicate(n.trials, {
  sigma <- runif(1, .03, .2)
  E <- 0
  Var <- sigma^2
  x <- rnorm(n.sample, sd = sigma)
  c(mean(x) - E, var(x) - Var)
})
deltas <- replicate(n.trials, {
  p <- runif(1)
  n <- sample(10:100, 1)
  E <- p
  Var <- p*(1-p) / n
  x <- rbinom(n.sample, n, p) / n
  c(mean(x) - E, var(x) - Var)
})
deltas <- replicate(n.trials, {

})


deltas <- replicate(n.trials, {
  l <- sample(1:100, 1)
  mu <- l
  x <- rpois(n.sample, l)
  mean(x) - mu
})

plot(density(deltas));abline(v=0)

stopImplicitCluster()
