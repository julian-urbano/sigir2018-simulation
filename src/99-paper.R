source("src/common.R")
source("src/io.R")

library(VineCopula)
library(simIReff)
library(corrplot)

measure.prettynames <- list(ap = "AP", ndcg20 = "nDCG@20", err20 = "ERR@20",
                            p10 = "P@10", p20 = "P@20", rr = "RR")
measure.col <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#a65628')
measure.pch <- c(1,2,3,4,6,22)
collection.prettynames <- list(web2010 = "Web 2010", web2011 = "Web 2011", web2012 = "Web 2012",
                               web2013 = "Web 2013", web2014 = "Web 2014")
collection.col <- c(web2010 = '#e41a1c', web2011 = '#377eb8', web2012 = '#4daf4a',
                     web2013 = '#984ea3', web2014 = '#ff7f00')
cont.col <- c("#6baed6", "#08519c", "#fd8d3c", "#a63603") #c('#bae4b3','#bdd7e7','#74c476','#6baed6')
disc.col <- c("#6a51a3", "#74c476", "#006d2c") #c('#bcbddc','#fdae6b','#e6550d')

# Examples =========================================================================================

my.dev.new("output/figs/cop.pdf", num = 6, ratio = 1)
x <- read.csv("data/web2010_ap.csv")[,c(69,16)]
x <- pobs(x, ties.method = "random")
plot(qnorm(x[,1]), qnorm(x[,2]), col = "grey", pch = 19, xlim = c(-3,3), ylim = c(-3,3),
       xlab = "System 1", ylab = "System 2", main = "AP", axes = FALSE)
my.axis(1, at=-3:3)
my.axis(2, at = -3:3)
box()
cop <- bicop(x, family_set = c("onepar","twopar"), selcrit = "loglik", presel = FALSE)
contour(cop, drawlabels = FALSE, col = collection.col[1], add = TRUE)
cop <- bicop(x, family_set = "gauss")
contour(cop, drawlabels = FALSE, col = collection.col[2], add = TRUE)
legend("topleft", cex = .8, bty = "n", col = collection.col[1:2], lwd = 1, legend = c("Clayton","Normal"))
my.dev.off()

my.dev.new("output/figs/cont.pdf", num = 6, ratio = 1)
x <- read.csv("data/web2010_ndcg20.csv")[,25]
x01 <- seq(1e-6, 1-1e-6, length.out = 100)
hist(x, xlim = 0:1, ylim = c(0,3.5), probability = TRUE,
     main = "Continuous", xlab = "nDCG@20", border = "grey")
box()
abline(h=0, col="grey")
effs <- effContFit(x)
sapply(seq_along(effs), function(e) lines(x01, deff(x01, effs[[e]]), col = cont.col[e]))
legend("topright", cex = .8, bty = "n", col = cont.col, lwd = 1, legend = c("N", "B", "NKS", "BKS"))
my.dev.off()

my.dev.new("output/figs/disc.pdf", num = 6, ratio = 1)
x <- read.csv("data/web2010_p10.csv")[,29]
x01 <- seq(0,1,.1)
plot(x01, tabulate(x*10+1,11) / length(x), type = "o", ylim = c(0,.3),
     main = "Discrete", xlab = "P@10", col = "grey", ylab = "Mass", pch = 19)
effs <- effDiscFit(x, x01)
sapply(seq_along(effs), function(e) lines(x01, deff(x01, effs[[e]]), col = disc.col[e], type = "o"))
legend("topright", cex = .8, bty = "n", col = disc.col, lwd = 1, legend = c("BB", "DKS", "DKS-2"))
my.dev.off()

# EVALUATION #######################################################################################

## Margins =========================================================================================

codeSelections <- function(dat, cont) {
  dat <- as.matrix(dat)

  dat <- gsub("t\\(([^\\)]+)\\)", "\\1", dat) # transformed
  dat <- gsub("dks(1)", "dks", dat, fixed = TRUE) # dks(1) = dks
  dat <- gsub("dks\\(.+\\)", "dks(m)", dat) # dks(m)

  effs <- c("norm", "beta", "nks", "bks")
  if(!cont) effs <- c("bbinom", "dks", "dks(m)")
  dat <- apply(dat, 2, function(d) table(factor(d, levels = effs)) / nrow(dat)*100)
  colnames(dat) <- c("LL", "AIC", "BIC")
  dat
}

readSelections <- function(trans, cont) {
  measures <- c("ap","ndcg20","err20")
  if(!cont) measures <- c("p10","p20","rr")
  dir <- ifelse(trans, "output/04-transform_select/", "output/03-margins_select/")

  dat.all <- NULL
  for(measure in measures) {
    dat <- NULL
    for(collection in .COLLECTIONS) {
      dat.f <- read.csv(paste0(dir, collection,"_",measure,".csv"), stringsAsFactors = FALSE)
      dat <- rbind(dat, dat.f)
    }
    dat <- codeSelections(dat, cont)
    dat.all <- cbind(dat.all, NA, dat)
  }
  dat.all <- dat.all[,-1]
  dat.all
}

w <- rep(c(1,1,1,.25),3)
s <- .2
at <- cumsum(mean(w)*s+w)-w/2
for(trans in c(TRUE, FALSE)) {
  for(cont in c(TRUE, FALSE)) {
    my.dev.new(file = paste0("output/figs/", ifelse(trans, "t", "u"),
                             "effdistSelect_", ifelse(cont, "cont", "disc"), ".pdf"),
               num = 4, ratio = .55)
    dat <- readSelections(trans, cont)

    if(cont) {
      barplot(dat, width = w, space = s, col = cont.col, xlab = "Selection criterion",
              ylab = "% cases", cex.names = .9, axes = FALSE, axisnames = FALSE)
      my.axis(1, at = at, colnames(dat), lwd=0)
      my.axis(2, c(0,25,50,75,100))

      title(measure.prettynames$ap, adj = .15)
      title(measure.prettynames$ndcg20, adj = .5)
      title(measure.prettynames$err20, adj = .91)
    }else{
      barplot(dat, width = w, space = s, col = disc.col, xlab = "Selection criterion",
              ylab = "% cases", cex.names = .9, axes = FALSE, axisnames = FALSE)
      my.axis(1, at = at, colnames(dat), lwd=0)
      my.axis(2, c(0,25,50,75,100))

      title(measure.prettynames$p10, adj = .14)
      title(measure.prettynames$p20, adj = .5)
      title(measure.prettynames$rr, adj = .85)
    }
    my.dev.off()
  }
}

my.dev.new(file = "output/figs/legend_cont.pdf", num = 4, ratio = .12, mar = c(-2.5,-2.5,-1.8,0))
plot(NA, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab = "", ylab = "")
legend("top", legend = c("Normal (N)", "Beta (B)", "Normal KS (NKS)", "Beta KS (BKS)"),
       fill = cont.col, bty = "n", ncol = 2 ,cex = .8)
my.dev.off()
my.dev.new(file = "output/figs/legend_disc.pdf", num = 4, ratio = .12, mar = c(-2.5,-2.5,-1.8,0))
plot(NA, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab = "", ylab = "")
legend("top", legend = c("Beta-Binomial (BB)", NA, "Discrete KS (DKS)", "DKS w/ mult. (DKS-h)"),
       fill = c(disc.col[1],NA,disc.col[-1]), bty = "n", ncol = 2, cex = .8)
my.dev.off()

## Copulas =========================================================================================

### Types of pair copulas --------------------------------------------------------------------------

families <- NULL
for(file in list.files("output/05-bicop_types/", full.names = TRUE)) {
  dat <- read.csv(file, stringsAsFactors = FALSE)
  families <- rbind(families, dat, stringsAsFactors = FALSE)
}
families <- BiCopName(families$logLik, FALSE) # logLik, AIC or BIC
families <- gsub("^(Survival|Rotated) ", "", families)
families <- gsub(" (90|180|270) degrees$", "", families)
families <- gsub("  ", " ", families)
tab <- as.data.frame(table(families))
sum(tab$Freq)
tab$Freq <- tab$Freq / sum(tab$Freq) * 100
tab[order(tab$Freq, decreasing = TRUE),]

### Goodness of fit --------------------------------------------------------------------------------

res <- expand.grid(collection = .COLLECTIONS, measure = .MEASURES,
                   logLik_rvine = 0, logLik_gaussian = 0,
                   AIC_rvine = 0, AIC_gaussian = 0,
                   BIC_rvine = 0, BIC_gaussian = 0,
                   n = 0, stringsAsFactors = FALSE)
for(i in 1:nrow(res)) {
  for(coptype in c("gaussian", "rvine")) {
    for(selcrit in c("logLik", "AIC", "BIC")) {
      cop <- load.object(file.path("scratch/06-cop_fit",
                                   paste0(res$collection[i], "_", res$measure[i]),
                                   paste0(selcrit, "_", coptype, "_full")))
      res[i, paste0(selcrit, "_", coptype)] <- get(selcrit)(cop)
      res[i, "n"] <- ncol(cop$matrix)
    }
  }
}

my.dev.new(file = "output/figs/logLik.pdf", num = 6, ratio = .9)
plot(res$logLik_gaussian, res$logLik_rvine, xlim=c(600,7700),ylim=c(600,7700),
     xlab = "Gaussian copula", ylab = "R-vine copula", main = "LL", col = collection.col)
abline(0:1)
my.dev.off()

my.dev.new(file = "output/figs/AIC.pdf", num = 6, ratio = .9)
plot(res$AIC_gaussian, res$AIC_rvine, xlim=c(-8500,-470),ylim=c(-8500,-470),
     xlab = "Gaussian copula", ylab = "R-vine copula", main = "AIC", col = collection.col)
abline(0:1)
my.dev.off()

my.dev.new(file = "output/figs/BIC.pdf", num = 6, ratio = .9)
plot(res$BIC_gaussian, res$BIC_rvine, xlim=c(-4200,3850),ylim=c(-4200,3850),
     xlab = "Gaussian copula", ylab = "R-vine copula", main = "BIC", col = collection.col)
abline(0:1)
my.dev.off()

my.dev.new(file = "output/figs/legend_collections.pdf",
           num = 2, ratio = .04, mar = c(-2.5,-2.5,-1.8,0))
plot(NA, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab = "", ylab = "")
legend("top", legend = collection.prettynames, lwd = 2, col = collection.col, bty = "n",
       ncol = 5, cex = .8)
my.dev.off()

### tau per vine level -----------------------------------------------------------------------------

for(fun in c("mean", "max")) {
  my.dev.new(paste0("output/figs/tau_", fun, ".pdf"), num = 4.5, ratio = .8)
  plot(NA, xlim=c(1,78), ylim=c(.05,1), log="xy", axes= FALSE,
       main = list(mean = "Mean", max = "Maximum")[[fun]],
       xlab = "Truncation Level", ylab = expression("|"*tau*"|"))
  my.axis(1, at = c(1,2,5,10,20,50))
  my.axis(2, at = c(.05,.1,.2,.5,1))
  box()

  FUN <- get(fun)
  for(collection in .COLLECTIONS) {
    taus <- NULL
    for(f in list.files("scratch/06-cop_fit/", "^AIC_rvine+",
                        recursive = TRUE, full.names = TRUE)) {
      if(grepl(collection, f)) {
        cat(f,"\n")
        cop <- load.object(f)

        tau <- get_all_ktaus(cop)
        tau <- sapply(tau, function(ttau) FUN(abs(unlist(ttau))))
        taus <- rbind(taus, data.frame(l = seq_along(tau), tau = tau))
      }
    }
    taus <- taus[complete.cases(taus),]
    points(taus, col = grDevices::adjustcolor(collection.col[collection], alpha.f = .1))
    s <- smooth.spline(taus$l, taus$tau)
    lines(s, lwd = 2, col=collection.col[collection])
  }
  my.dev.off()
}

## Simulation ======================================================================================

### Correlation matrix -----------------------------------------------------------------------------

collection <- "web2010"
measure <- "ndcg20"

dat <- read.csv(paste0("data/", collection,"_",measure,".csv"))
dup.i <- indexOfDuplicates(dat)
if(length(dup.i) > 0)
  dat <- dat[, -dup.i]
effs <- lapply(colnames(dat), function(sysname)
  load.object(file.path("scratch/03-margins_select/", paste0(collection, "_", measure), "AIC",
                        sysname)))

plotcorr <- function(file , M, title) {
  my.dev.new(file=file, num=5, ratio = .97, mgp=-c(1,.7,0)/2, mar=-c(2.2,2,-.09,-1.95))
  corrplot(M, method = "color", title = title, tl.pos = "n", mar=c(0,0,1.6,0), cl.cex = .8)
  box()
  my.dev.off()
}

sim <- function(n, cop, effs) {
  r <- rvinecop(n, cop)
  sapply(1:ncol(r), function(i) qeff(r[,i], effs[[i]]))
}

set.seed(1)
copg <- load.object(file.path("scratch/06-cop_fit/",
                              paste0(collection, "_", measure), "logLik_gaussian_full"))
copr <- load.object(file.path("scratch/06-cop_fit/",
                              paste0(collection, "_", measure), "logLik_rvine_full"))
copr2 <- vinecop(pobs(dat, ties.method = "random"), family_set = c("onepar", "twopar"),
                 selcrit = "loglik", cores = max(1, detectCores()-1), presel = FALSE, trunc_lvl = 2)

simg <- sim(500, copg, effs)
simr <- sim(500, copr, effs)
simr2 <- sim(500, copr2, effs)

M.simg <- cor(simg, method = "s")
M.simr <- cor(simr, method = "s")
M.simr2 <- cor(simr2, method = "s")
M.dat <- cor(dat, method = "s")

i <- corrMatOrder(M.dat, order = "hclust")
plotcorr("output/figs/cor_dat.pdf", M.dat[i,i], title = "Original data")
plotcorr("output/figs/cor_g.pdf", M.simg[i,i], title = "Gaussian copula")
plotcorr("output/figs/cor_r.pdf", M.simr[i,i], title = "R-vine copula (no truncation)")
plotcorr("output/figs/cor_r2.pdf", M.simr2[i,i], title = "R-vine copula (2 levels)")

# Copula Precision ---------------------------------------------------------------------------------

readData <- function(measure, w) {
  sapply(.COLLECTIONS, function(collection) {
    r <- read.csv(file.path("output/07-cop_sim/", paste0(collection, "_", measure, "_", w, ".csv")))
    effs <- lapply(colnames(r), function(rr)
      load.object(file.path("scratch/03-margins_select", paste0(collection, "_", measure),
                            "AIC", rr)))
    truth <- sapply(effs, function(eff) eff[[w]])
    truth <- matrix(truth, byrow = TRUE, ncol = ncol(r), nrow = nrow(r))
    unlist(r - truth)
  })
}

for(w in c("mean", "var")) {
  d_ap <- density(unlist(readData("ap", w)))
  d_ndcg20 <- density(unlist(readData("ndcg20", w)))
  d_err20 <- density(unlist(readData("err20", w)))
  d_p10 <- density(unlist(readData("p10", w)))
  d_p20 <- density(unlist(readData("p20", w)))
  d_rr <- density(unlist(readData("rr", w)))

  my.dev.new(paste0("output/figs/devs_cop_", w, ".pdf"), num = 4, ratio = .7,
             mar = c(0,-2,0,.7))
  plot(NA, xlim = c(-1,1)*list(mean = .025, var = .006)[[w]],
       # max(abs(c(d_ap$x, d_ndcg20$x, d_err20$x, d_p10$x, d_p20$x, d_rr$x))),
       ylim = c(0, max(d_ap$y, d_ndcg20$y, d_err20$y, d_p10$y, d_p20$y, d_rr$y)),
       axes = FALSE, ylab="", xlab = list(mean = expression("Deviation from " * mu),
                                          var = expression("Deviation from " * sigma^2))[[w]],
       main = "Copula")
  abline(v=0)
  my.axis(1, -1:1)
  if(w == "mean") my.axis(1, at = seq(-.02,.02,.01))
  if(w == "var")  my.axis(1, at = seq(-.006,.006,.002))

  lines(d_ap, col = measure.col[1])
  lines(d_ndcg20, col = measure.col[2])
  lines(d_err20, col = measure.col[3])
  lines(d_p10, col = measure.col[4])
  lines(d_p20, col = measure.col[5])
  lines(d_rr, col = measure.col[6])

  if(w=="mean")
    legend("topleft", legend = measure.prettynames, lwd = 2, col = measure.col, bty = "n", cex= .75)

  my.dev.off()
}

# R Precision --------------------------------------------------------------------------------------

set.seed(1)
trials <- 100000
d_beta <- replicate(trials, {
  shape1 <- runif(1, 1, 20)
  shape2 <- runif(1, 1, 20)

  x <- rbeta(1000, shape1, shape2)
  E <- shape1 / (shape1 + shape2)
  Var <- shape1 * shape2 / (shape1 + shape2)^2 / (shape1 + shape2 + 1)
  c(mean = mean(x)-E, var = var(x)-Var)
})
d_norm <- replicate(trials, {
  sigma <- runif(1, .15, .2)

  x <- rnorm(1000, sd = sigma)
  E <- 0
  Var <- sigma^2
  c(mean = mean(x)-E, var = var(x)-Var)
})
d_bbinom <- replicate(trials, {
  n <- 10
  shape1 <- runif(1, 1, 6)
  shape2 <- runif(1, 1, 6)

  x <- extraDistr::rbbinom(1000, n, shape1, shape2) / n
  E <- shape1 / (shape1 + shape2)
  Var <- shape1 * shape2 * (shape1 + shape2 + n) / (shape1 + shape2)^2 / (shape1 + shape2 + 1) / n
  c(mean = mean(x)-E, var = var(x)-Var)
})

for(w in c("mean", "var")) {
  my.dev.new(paste0("output/figs/devs_R_", w, ".pdf"), num = 4, ratio = .7,
             mar = c(0,-2,0,.7))
  plot(NA, xlim = c(-1,1)*list(mean = .025, var = .006)[[w]],
       ylim = c(0, list(mean=130,var=950)[[w]]),
       axes = FALSE, ylab="", xlab = list(mean = expression("Deviation from " * mu),
                                          var = expression("Deviation from " * sigma^2))[[w]],
       main = "R")
  abline(v=0)
  my.axis(1, -1:1)
  if(w == "mean") my.axis(1, at = seq(-.02,.02,.01))
  if(w == "var")  my.axis(1, at = seq(-.006,.006,.002))

  lines(density(d_beta[w,]), col = measure.col[1])#, lty=2)
  lines(density(d_norm[w,]), col = measure.col[2])#, lty=2)
  lines(density(d_bbinom[w,]), col = measure.col[3])#, lty=2)

  if(w == "mean")
    legend("topleft", legend = c("Beta","Normal","Beta-\nBinomial"), lwd = 2,
           col = measure.col[1:3], bty = "n", cex = .75)

  my.dev.off()
}

# SAMPLE APPLICATIONS ##############################################################################

# Statistical Power and Topic Set Size =============================================================

lim <- c(.005,.2)

dat <- read.csv("scratch/10-power/res.csv")
dat$sigma.true <- 1
a <- aggregate(cbind(sigma.true, sigma.150, sigma.hat, ntopics)~trial, dat, mean)
my.dev.new(file = paste0("output/figs/power_resample.pdf"), num = 4.5, ratio = .8)
plot(a$sigma.150, a$sigma.hat, xlim = lim, ylim = lim,
     xlab = expression(sigma[150]), ylab = expression(sigma[B]), main = "O) Resampling")
abline(lm(sigma.hat~sigma.150,a),col=2,lty=2)
abline(0:1,col=2)
my.dev.off()
mean(a$sigma.hat / a$sigma.150)
mean(a$ntopics)

dat <- read.csv("scratch/10-power/cop_res.csv")
a <- aggregate(cbind(sigma.true, sigma.150, sigma.hat, ntopics)~trial, dat, mean)
my.dev.new(file = paste0("output/figs/power_cop_res1.pdf"), num = 4.5, ratio = .8)
plot(a$sigma.150, a$sigma.hat, xlim = lim, ylim = lim,
     xlab = expression(sigma[150]), ylab = expression(sigma[B]), main = "A1) Copula + Resampling")
abline(lm(sigma.hat~sigma.150,a),col=2,lty=2)
abline(0:1,col=2)
my.dev.off()

my.dev.new(file = paste0("output/figs/power_cop_res2.pdf"), num = 4.5, ratio = .8)
plot(a$sigma.true, a$sigma.hat, xlim = lim, ylim = lim,
     xlab = expression(sigma), ylab = expression(sigma[B]), main = "A1) Copula + Resampling")
abline(lm(sigma.hat~sigma.true,a),col=2,lty=2)
abline(0:1,col=2)
my.dev.off()
mean(a$sigma.hat / a$sigma.true)
mean(a$ntopics)

dat <- read.csv("scratch/10-power/cop.csv")
a <- aggregate(cbind(sigma.true, sigma.150, sigma.hat, ntopics)~trial, dat, mean)
my.dev.new(file = paste0("output/figs/power_cop.pdf"), num = 4.5, ratio = .8)
plot(a$sigma.true, a$sigma.hat, xlim = lim, ylim = lim,
     xlab = expression(sigma), ylab = expression(sigma[B]), main = "A2) Copula")
abline(lm(sigma.hat~sigma.true,a),col=2,lty=2)
abline(0:1,col=2)
my.dev.off()
mean(a$sigma.hat / a$sigma.true)
mean(a$ntopics)

# Hypothesis Testing and Type I Errors =============================================================

# helper function to code the conflicts found in the experiment's result
codeConflicts <- function(r) {
  r[r %in% c("<<d")] <- "major"
  r[r %in% c("<>d", "><d")] <- "minor"
  r[r %in% c("<<e")] <- "success"
  r[r %in% c("<>e","><e")] <- "lack"
  r[r %in% c(">>e",">>d")] <- "inconclusive"
  as.list(table(r))
}

for(measure in c("ap", "p10")) {
  # Split-half
  r <- read.csv(paste0("scratch/11-type_I/", measure, "_splithalf.csv"))$res
  # every other row (a pair (i,j) appears also as (j,i) right after)
  r <- as.character(r[seq(2,length(r), by = 2)])
  r <- codeConflicts(r)
  
  signif <- r$success*2 + r$major*2 + r$minor + r$lack
  cat("splithalf", measure, (r$major*2+r$minor)/signif*100, "\n")
  
  # cop_splithalf
  r <- read.csv(paste0("scratch/11-type_I/", measure, "_cop_splithalf.csv"))$res
  # every other row (a pair (i,j) appears also as (j,i) right after)
  r <- as.character(r[seq(2,length(r), by = 2)])
  r <- codeConflicts(r)
  
  signif <- r$success*2 + r$major*2 + r$minor + r$lack
  cat("cop_splithalf", measure, (r$major*2+r$minor)/signif*100, "\n")
  
  # cop_same
  r <- read.csv(paste0("scratch/11-type_I/", measure, "_cop_same.csv"))$p
  n <- sum(!is.na(r))
  cat("cop_same .05", measure, sum(r<.05, na.rm = TRUE) / n*100, "\n")
  cat("cop_same .01", measure, sum(r<.01, na.rm = TRUE) / n*100, "\n")
  
  # cop_trans
  r <- read.csv(paste0("scratch/11-type_I/", measure, "_cop_trans.csv"))$p
  n <- sum(!is.na(r))
  cat("cop_trans .05", measure, sum(r<.05, na.rm = TRUE) / n*100, "\n")
  cat("cop_trans .01", measure, sum(r<.01, na.rm = TRUE) / n*100, "\n")
}
