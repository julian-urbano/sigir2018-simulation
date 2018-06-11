source("src/common.R")
library(rvinecopulib)
library(VineCopula)
library(doParallel)
stopImplicitCluster()
registerDoParallel(cores = max(1, detectCores()-1))

# Helper functions to run iterative sampling -------------------------------------------------------

# use sampleTopic.FUN to iteratively sample topics until 80% power is achieved to detect an effect
# of delta.target at alpha=.05. A minimum of min.ntopics is required.
iterativeSampling <- function(delta.target, min.ntopics, sampleTopic.FUN, ...) {
  x <- sampleTopic.FUN(n = 1000, ...) # pre-sample many between-system deltas for efficiency
  ntopics <- min.ntopics-1
  mu.hat <- NA
  sigma.hat <- NA

  delta.detectable <- Inf # detectable delta
  while(delta.detectable > delta.target) {
    ntopics <- ntopics + 1
    mu.hat <- mean(x[1:ntopics])
    sigma.hat <- sd(x[1:ntopics])

    if(!is.na(sigma.hat) && sigma.hat > 0) # beware of cases with no variability
      delta.detectable <- power.t.test(n = ntopics, sd = sigma.hat, power = .8, sig.level = .05,
                        type = "one", alternative = "two")$delta
  }
  return(c(ntopics = ntopics, mu.hat = mu.hat, sigma.hat = sigma.hat,
           delta.detectable = delta.detectable))
}

# to resample n topics from D
sampleTopic.resample <- function(n, D) {
  sample(D, size = n, replace = TRUE)
}
# to sample n topics from the copula and margins
sampleTopic.copula <- function(n, effA, effB, cop) {
  u <- rbicop(n, cop)
  qeff(u[,1], effA) - qeff(u[,2], effB)
}

# Actual experiment --------------------------------------------------------------------------------

out.path <- "scratch/10-power"
dir.create(out.path, recursive = TRUE, showWarnings = FALSE)

# Read Robust 2004 (no Description-only runs, topics 301-450) and sort by mean score
dat <- read.csv("data/robust2004(noD)_ap.csv")
dat <- dat[, order(colMeans(dat), decreasing = TRUE)]

# 'res' is for resampling (design O in the paper), 'cop_res' to resample from the copula
# (design A1), and 'cop' to sample from the copula (design A2)
for(sample.FUN.name in c("res", "cop_res", "cop")) {
  trials1 <- 200 # number of system pairs to use
  trials2 <- 1000 # number of trials per system pair
  res <- foreach(trial = 1:trials1, .combine = rbind) %dopar% {
    cat("\n", sample.FUN.name, trial);flush.console()
    
    set.seed(trial)
    # Select systems at random and compute true mean and sigma
    dat.quartile <- as.integer(cut(seq_along(dat), 4))
    A.name <- sample(names(dat)[dat.quartile %in% 1:3], size = 1) # top 3 quartiles
    B.name <- sample(names(dat)[dat.quartile == 2], size = 1) # 2nd quartile
    A <- dat[,A.name]
    B <- dat[,B.name]
    D <- A-B
    
    if(A.name != B.name) {
      if(sample.FUN.name == "res") {
        mu.true <- mean(D)
        sigma.150 <- sd(D)
        sigma.true <- NA # unknown
        sigma.delta <- sigma.150 # to compute target delta
      }else if(sample.FUN.name == "cop_res"){
        effA <- effContFitAndSelect(A)
        effB <- effContFitAndSelect(B)
        cop <- bicop(cbind(pobs(A, ties.method = "random"), pobs(B, ties.method = "random")),
                     family_set = c("onepar", "twopar")) # fit bivariate coputa
        mu.true <- effA$mean - effB$mean
        # for simplicity, approximate true sigma via monte carlo
        sigma.true <- mean(replicate(1000, sd(sampleTopic.copula(200, effA, effB, cop))))
        D <- sampleTopic.copula(150, effA, effB, cop)
        sigma.150 <- sd(D)
        sigma.delta <- sigma.true # to compute target delta
      } else if(sample.FUN.name == "cop") {
        effA <- effContFitAndSelect(A)
        effB <- effContFitAndSelect(B)
        cop <- bicop(cbind(pobs(A, ties.method = "random"), pobs(B, ties.method = "random")),
                     family_set = c("onepar", "twopar")) # fit bivariate copula
        mu.true <- effA$mean - effB$mean
        # for simplicity, approximate true sigma via monte carlo
        sigma.true <- mean(replicate(1000, sd(sampleTopic.copula(200, effA, effB, cop))))
        sigma.150 <- sd(sampleTopic.copula(150, effA, effB, cop))
        sigma.delta <- sigma.true # to compute target delta
      }
      
      if(sigma.delta > 0) {
        # this is the desired delta to detect after iterative sampling
        delta.target <- abs(power.t.test(n = 100, sd = sigma.delta, power = .8, sig.level = .05,
                                         type = "one", alternative = "two")$delta)
        # sample topics until we reach enough power to detect delta
        if(sample.FUN.name == "res")
          dat.trial <- replicate(trials2,
                                 iterativeSampling(delta.target, 40, sampleTopic.resample,
                                                   D = D))
        else if(sample.FUN.name == "cop")
          dat.trial <- replicate(trials2,
                                 iterativeSampling(delta.target, 40, sampleTopic.copula,
                                                   effA = effA, effB = effB, cop = cop))
        else if(sample.FUN.name == "cop_res")
          dat.trial <- replicate(trials2,
                                 iterativeSampling(delta.target, 40, sampleTopic.resample,
                                                   D = D))
        else stop("sample.FUN.name?")
        
        cbind(trial = trial, mu.true = mu.true, sigma.true = sigma.true, sigma.150 = sigma.150,
              delta.target = delta.target, t(dat.trial)) # append all data
      }
    }
  }
  write.csv(file = file.path(out.path, paste0(sample.FUN.name, ".csv")), row.names = FALSE, res)
}

stopImplicitCluster()
