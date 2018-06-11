source("src/common.R")
library(VineCopula)
library(rvinecopulib)
library(doParallel)
stopImplicitCluster()
registerDoParallel(cores = max(1, detectCores()-1))

# helper function to code the result of a split-half comparison
replicaResult <- function(d1, d2) {
  mu1 <- mean(d1)
  mu2 <- mean(d2)
  p1 <- t.test(d1)$p.value
  p2 <- t.test(d2)$p.value
  
  l1 <- ifelse(p1<.05, "<",">")
  l2 <- ifelse(p2<.05, "<",">")
  l3 <- ifelse(sign(mu1)==sign(mu2), "e","d")
  return(paste0(l1,l2,l3))
}

datAP <- read.csv("data/robust2004_ap.csv") # all original AP scores, to subset systems later on
out.path <- "scratch/11-type_I"
dir.create(out.path, showWarnings = FALSE, recursive = TRUE)

for(measure in c("ap", "p10")) {
  # read original data, and subset TREC 7-8 topics, only top 75% of systems over all topics
  dat <- read.csv(paste0("data/robust2004_", measure, ".csv"))
  dat <- dat[dat$topic %in% 351:450, indexOfTopN(datAP[-1], p = .75)+1]
  n <- ncol(dat)
  nn <- n*(n-1)/2
  
  # Split-half resampling from the original data (design O in the paper) ===========================
  
  trials <- 1000
  # iterate over trials
  allres <- foreach(trial = 1:trials, .combine = rbind) %dopar% {
    cat("\n",trial);flush.console()
    set.seed(trial)
    ti <- sample(1:nrow(dat), size = 50) # topic ids of first half
    
    res <- data.frame(i = rep(0, nn*2), j = 0, trial = 0, res = "", stringsAsFactors = FALSE)
    res.i <- 0
    for(i in 1:(n-1)) { # for every system index i
      si <- colnames(dat)[i]
      for(j in (i+1):n){ # for every other system index j
        sj <- colnames(dat)[j]
        d <- dat[,i]-dat[,j]
        
        res.i <- res.i +1
        res[res.i,] <- list(i = i, j = j, trial = trial, res = replicaResult(d[ti],d[-ti]))
        res.i <- res.i +1
        res[res.i,] <- list(i = i, j = j, trial = trial, res = replicaResult(d[-ti],d[ti]))
      }
    }
    res # return to be collected with the other trials
  }
  write.csv(file = file.path(out.path, paste0(measure, "_splithalf.csv")), row.names = FALSE,
            allres)
  
  # Split-half resampling from the copula (design A1) ==============================================
  
  trials <- 1000
  allres <- foreach(i = 1:(n-1), .combine = rbind) %dopar% { # iterate over systems
    library(simIReff)
    library(VineCopula)
    library(rvinecopulib)
    cat("\n",i);flush.console()
    set.seed(i)
    
    # margin for system i
    if(measure == "ap")
      effi <- effContFitAndSelect(dat[,i])
    else
      effi <- effDiscFitAndSelect(dat[,i], support(measure))
    
    res <- data.frame(i = rep(0, (n-i)*trials*2), j = 0, trial = 0, res = "",
                      stringsAsFactors = FALSE)
    res.i <- 0
    for(j in (i+1):n) { # for every other system j
      if(abs(cor(dat[,i], dat[,j], method = "k")) < .99) { # only for systems with variability
        # margin for system j
        if(measure == "ap")
          effj <- effContFitAndSelect(dat[,j])
        else
          effj <- effDiscFitAndSelect(dat[,j], support(measure))
        cop <- bicop(pobs(dat[,c(i,j)], ties.method = "random"),
                     family_set = c("onepar","twopar"), presel = FALSE)
        
        for(trial in 1:trials) {
          u <- rbicop(100, cop)
          d <- qeff(u[,1], effi) - qeff(u[,2], effj)
          
          res.i <- res.i +1
          res[res.i,] <- list(i = i, j = j, trial = trial, res = replicaResult(d[1:50],d[51:100]))
          res.i <- res.i +1
          res[res.i,] <- list(i = i, j = j, trial = trial, res = replicaResult(d[51:100],d[1:50]))
        }
      }
    }
    res # return to be collected with the other trials
  }
  write.csv(file = file.path(out.path, paste0(measure, "_cop_splithalf.csv")), row.names = FALSE,
            allres)
  
  # Simulation from the copula with same margins (design A2) =======================================
  
  trials <- 10000
  allres <- foreach(i = 1:n, .combine = rbind) %dopar% { # iterate over systems
    library(simIReff)
    library(VineCopula)
    library(rvinecopulib)
    cat("\n",i);flush.console()
    set.seed(i)
    
    # margin for system i
    if(measure == "ap")
      effi <- effContFitAndSelect(dat[,i])
    else
      effi <- effDiscFitAndSelect(dat[,i], support(measure))
    
    # for every other candidate system j, fit a bivariate copula with i
    copsj <- lapply(1:n, function(j) {
      if(abs(cor(dat[,i], dat[,j], method = "k")) < .99)
        bicop(pobs(dat[,c(i,j)], ties.method = "random"),
              family_set = c("onepar","twopar"), presel = FALSE)
      else
        NULL
    })
    
    res <- data.frame(i = rep(i, trials), j = 0, trial = 0, p = NA, stringsAsFactors = FALSE)
    for(trial in 1:trials) {
      j <- sample(seq_along(copsj), 1) # sample one system and get its bivariate copula with i
      copj <- copsj[[j]]
      
      if(!is.null(copj)) {
        u <- rbicop(50, copj)
        
        xi <- qeff(u[,1], effi)
        xj <- qeff(u[,2], effi) # same margin as i
        
        p <- t.test(xi-xj)$p.value
        res[trial,] <- list(i = i, j = j, trial = trial, p = p)
      }else{
        res[trial,] <- list(i = i, j = j, trial = trial, p = NA)
      }
    }
    res
  }
  write.csv(file = file.path(out.path, paste0(measure, "_cop_same.csv")), row.names = FALSE, allres)
  
  # Simulation from the copula with transformed margins (design A3) ================================
  
  trials <- 10000
  mu <- colMeans(dat) # means observed on original data
  
  allres <- foreach(i = 1:n, .combine = rbind) %dopar% { # iterate over systems
    library(simIReff)
    library(VineCopula)
    library(rvinecopulib)
    cat("\n",i);flush.console()
    set.seed(i)
    
    # margin for system i
    if(measure == "ap")
      effi <- effContFitAndSelect(dat[,i])
    else
      effi <- effDiscFitAndSelect(dat[,i], support(measure))
    
    # candidate systems to compare with: closest 10 w.r.t. mean score
    js <- order(abs(mu-mean(dat[,i])))[2:11] # remove 1st because it's i itself
    
    # for each candidate system j, fit margin and transform to i's mean
    effsj <- lapply(js, function(j) {
      if(measure == "ap")
        effs <- effContFit(dat[,j])
      else
        effs <- effDiscFit(dat[,j], support(measure))
      effs <- effTransformAll(effs, means = effi$mean)
      effs <- effs[!is.na(effs)] # reject those without successful transformation
      if(length(effs)>0)
        effSelect(effs)
      else
        NULL
    })
    
    # for eah successful candidate j, fit a bivariate copula with i
    copsj <- lapply(js, function(j) {
      if(abs(cor(dat[,i], dat[,j], method = "k")) < .99)
        bicop(pobs(dat[,c(i,j)], ties.method = "random"),
              family_set = c("onepar","twopar"), presel = FALSE)
      else
        NULL
    })    
    
    res <- data.frame(i = rep(i, trials), j = 0, trial = 0, p = NA, stringsAsFactors = FALSE)
    for(trial in 1:trials) {
      j_ <- sample(seq_along(js), 1) # index of j within js
      j <- js[j_] # actual index j
      
      # select transformed margin and copula
      effj <- effsj[[j_]]
      copj <- copsj[[j_]]
      
      if(!is.null(effj) & !is.null(copj)) {
        u <- rbicop(50, copj)
        
        xi <- qeff(u[,1], effi)
        xj <- qeff(u[,2], effj)
        
        p <- t.test(xi-xj)$p.value
        res[trial,] <- list(i = i, j = j, trial = trial, p = p)
      }else{
        res[trial,] <- list(i = i, j = j, trial = trial, p = NA)
      }
    }
    res
  }
  write.csv(file = file.path(out.path, paste0(measure, "_cop_trans.csv")), row.names = FALSE,
            allres)
}

stopImplicitCluster()
