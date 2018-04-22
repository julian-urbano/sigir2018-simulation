.COLLECTIONS <- c("web2010", "web2011", "web2012", "web2013", "web2014")
.MEASURES <- c("ap", "ndcg20", "err20", "p10", "p20", "rr")

# Set checkpoint and load simIReff #################################################################

checkpoint::checkpoint("2018-04-20")
if(!"simIReff" %in% installed.packages()[,1])
	devtools::install_github(repo = "julian-urbano/simIReff", ref = "7619a02cb9a4d5079e6afeb801a7b199d4c37109")

# Utils ############################################################################################

indexOfTopN <- function(X, p = .75) {
  mu <- colMeans(X)
  n <- ceiling(ncol(X)*p)
  i <- order(mu, decreasing = TRUE)
  return(i[1:n])
}

indexOfDuplicates <- function(dat) {
  toremove <- integer(0)
  for(i in seq(ncol(dat)-1)) {
    X <- dat[,i]
    for(j in seq(i+1,ncol(dat))) {
      Y <- dat[,j]
      if(all(X==Y))
        toremove <- c(toremove, j)
    }
  }
  return(toremove)
}

# Plots ############################################################################################

my.dev.width <- function(num=1){
  return(12 / num)
}
my.dev.par <- function(mar = 0, mgp = 0, ...){
  par(mar = c(3.1,2.9,2,0.2) + mar, mgp = c(1.9,.7,0) + mgp, ...)
}
my.dev.abline <- function(col="darkgrey", lwd=1, lty=2, ...){
  abline(col=col, lwd=lwd, lty=lty, ...)
}
my.dev.set.pdf <- function() {
  .GlobalEnv$my.dev.new <- function(file, num, ratio=.82, ...){
    width <- my.dev.width(num)
    height <- width*ratio
    pdf(file=file, width=width, height=height)
    my.dev.par(...)
  }
  .GlobalEnv$my.dev.off <- function(...) { off <- capture.output(dev.off(...)) }
}
my.dev.set.win <- function() {
  .GlobalEnv$my.dev.new <- function(file, num, ratio=.82, ...){
    width <- my.dev.width(num)
    height <- width*ratio
    #dev.new(width=width, height=height)
    my.dev.par(...)
  }
  .GlobalEnv$my.dev.off <- function(){}
}
my.dev.set.pdf()
#my.dev.set.win()
