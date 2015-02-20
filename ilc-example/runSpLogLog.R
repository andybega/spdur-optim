# Optimizing splolog in rcpp

if(Sys.info()['user']=='janus829' | Sys.info()['user']=='s7m'){
  setwd("~/Research/WardProjects/spduration-optim/ilc-example")}

loadPkg=function(toLoad){
  for(lib in toLoad){
  if(! lib %in% installed.packages()[,1])
    { install.packages(lib, repos='http://cran.rstudio.com/') }
  suppressMessages( library(lib, character.only=TRUE) ) }
}

# Load libraries
packs=c("corpcor", 'lineprof', 'microbenchmark')
loadPkg(packs)

# Load core estimation functions
core_funcs <- list.files("../core-functions", pattern="*.R")
for(func in core_funcs){source(paste0('../core-functions/', func))}

# Load example data
load("irc-data-mod.rda")

# Theme 1 with all data
model1 <- function() {
  spdur(duration ~ 1 + log10(i.matl.conf.DIStGOV.l1+1) + 
          log10(i.matl.coop.GOVtGOV.l1+1),
        atrisk ~ 1 + ldr.irregular + ldr.foreign + log10(mths.in.power+1),
        data=irc.data, silent=TRUE)
}

# Run time; use times>1 to eliminate overhead variance
# microbenchmark(model1(), times=1, unit="s")

# # Line by line profiling; might give different results each run
# m1_l1 <- lineprof({
#   spdur(duration ~ 1 + log10(i.matl.conf.DIStGOV.l1+1) + 
#           log10(i.matl.coop.GOVtGOV.l1+1),
#         atrisk ~ 1 + ldr.irregular + ldr.foreign + log10(mths.in.power+1),
#         data=irc.data, silent=TRUE)
# }, interval=0.1)

# shine(m1_l1)


# Look only at likelihood calculation
# (Y, X, Z saved from debugging spdur)
load("XYZ-matrices.rda")

# weib_lnl
# spweib_lnl

set.seed(6886)
theta <- rnorm(ncol(X) + ncol(Z) + 1)
# spweib_l1 <- lineprof(spweib_lnl(theta, Y, X, Z), torture=T)

# Testing out stuff
library(Rcpp)
library(RcppArmadillo)
library(rbenchmark)

sourceCpp("../core-functions/loglog_lnlC.cpp")

test=loglog_lnlCpp(theta, Y, X)
# sum(exp(-X %*% theta[1:ncol(X)]) == test)/length(test)


#' Regular Log-logistic likelihood
theta=theta; y=Y; X=X
# loglog_lnl <- function(theta, y, X) {
  beta <- theta[1:ncol(X)]
  g <- theta[ncol(X) + 1]
  d <- y[,1]
  ti <- y[,2]
  t0 <- y[,4]
  ly <- y[,3]
  lambda <- exp(-X%*%beta)
  alpha <- exp(-g)
  # ln.ft <- log(alpha) + alpha*log(lambda) + (alpha-1)*log(ti) -2 * log(1+(lambda*ti)^alpha)
  ln.ft <- (lambda*ti)^alpha
 
head(cbind(test, ln.ft))
sum(test==ln.ft)/length(test)

#   ln.st  <- -log((1+(lambda*ti)^alpha))
#   ln.st0 <- -log((1+(lambda*t0)^alpha))
#   cens <- ifelse((d==1) & (ly==0) | (d==0) , ln.st - ln.st0 , 0)
#   nocens <- ifelse((d==1) & (ly==1) , ln.ft - ln.st0 , 0)
#   logl <- sum(cens+nocens)
#   return(-logl)
# }