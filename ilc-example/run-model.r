#
#   Example code from ILC 2014 papers
#

setwd("~/Work/spduration-optim/ilc-example")

library(corpcor)
library(lineprof)
library(microbenchmark)

# Load core estimation functions
core_funcs <- list.files("../core-functions", pattern="*.R")
sapply(paste0("../core-functions/", core_funcs), source, .GlobalEnv)

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
microbenchmark(model1(), times=1, unit="s")

# Line by line profiling; might give different results each run
m1_l1 <- lineprof({
  spdur(duration ~ 1 + log10(i.matl.conf.DIStGOV.l1+1) + 
          log10(i.matl.coop.GOVtGOV.l1+1),
        atrisk ~ 1 + ldr.irregular + ldr.foreign + log10(mths.in.power+1),
        data=irc.data, silent=TRUE)
}, interval=0.1)

shine(m1_l1)


# Look only at likelihood calculation
# (Y, X, Z saved from debugging spdur)
load("XYZ-matrices.rda")

# weib_lnl
# spweib_lnl

theta <- rnorm(ncol(X) + ncol(Z) + 1)
spweib_l1 <- lineprof(spweib_lnl(theta, Y, X, Z), torture=T)

shine(spweib_l1)
# The slowest portion of the code are the two ifelse cens/nocens calculations

# Here's the overhead part of the code
# Parameters
rx <- ncol(X)
rz <- ncol(Z)
beta <- theta[1:rx]
gamma <- theta[(rx+1):(rx+rz)]
p <- theta[rx+rz+1]
lambda <- exp(-X%*%beta)
alpha <- exp(-p)
  
# Dependent variables
rc <- Y[,1]    # right-censored spell (1 if uncensored)
ti <- Y[,2]
t0 <- Y[,4]    # duration at previous period
lo <- Y[,3]    # last observation in spell
  
# Logistic for population split (1 = at risk)
pr1 <- plogis(Z%*%gamma)
pr0 <- plogis(Z%*%gamma, lower.tail=FALSE)
  
# Weibull 
#ln.ft <- log(alpha * lambda) + (alpha-1) * log(lambda * ti) - (lambda * ti)^alpha
ln.ft <- log(alpha) + (alpha) * log(lambda) + (alpha - 1) * log(ti) - (lambda * ti)^alpha
st    <- exp( -(lambda * ti)^alpha )    # S(t) for current period
st0   <- exp( -(lambda * t0)^alpha )    # S(t) for previous period

cens1 <- function() {
  # Likelihood contributions
  cens   <- ifelse((lo==0) | (rc==0), log(pr0 + (pr1 * st)) - log(pr0 + (pr1 * st0)), 0)
  nocens <- ifelse((rc==1) & (lo==1), log(pr1) + ln.ft      - log(pr0 + (pr1 * st0)), 0)

  logl <- sum(cens+nocens)
  return(-logl)
}
 
cens2 <- function() {
  cens_ind <- (lo==0) | (rc==0)
  logl <- vector("numeric", length=length(ti))
  
  # Likelihood contributions
  logl[cens_ind]  <- (log(pr0 + (pr1 * st)) - log(pr0 + (pr1 * st0)))[cens_ind]
  logl[!cens_ind] <- (log(pr1) + ln.ft      - log(pr0 + (pr1 * st0)))[!cens_ind]

  return(-sum(logl))
}

cens3 <- function() {
  logl <- vector("numeric", length=length(ti))
  
  # Likelihood contributions
  logl[(lo==0) | (rc==0)]  <- (log(pr0 + (pr1 * st)) - log(pr0 + (pr1 * st0)))[(lo==0) | (rc==0)]
  logl[!((lo==0) | (rc==0))] <- (log(pr1) + ln.ft      - log(pr0 + (pr1 * st0)))[!((lo==0) | (rc==0))]
  
  -sum(logl)
}

cens4 <- function() {
  # Likelihood contributions
  logl <- ifelse((lo==0) | (rc==0),
                 log(pr0 + (pr1 * st)) - log(pr0 + (pr1 * st0)),
                 log(pr1) + ln.ft      - log(pr0 + (pr1 * st0)))
  
  logl <- sum(logl)
  -logl
}



microbenchmark(cens1, cens2, cens3, cens4, times=10^5)

cens5 <- function(cens, pr0, pr1, st, st0, ln.ft) {
  # just break up the ifelse so we can see whether it or the calculations take 
  # so long
  # Likelihood contributions
  cens   <- log(pr0 + (pr1 * st)) - log(pr0 + (pr1 * st0))
  nocens <- log(pr1) + ln.ft      - log(pr0 + (pr1 * st0))
  
  logl <- ifelse((lo==0) | (rc==0), cens, nocens)
  -sum(logl)
}

c5_l <- lineprof(cens5, torture=T)


library(Rcpp)

cppFunction('double spweib_lnlC(LogicalVector c, NumericVector pr0, NumericVector pr1,
  NumericVector st, NumericVector st0, NumericVector ln_ft) {
  int n = c.size();
  double lnl = 0;
  
  for(int i = 0; i < n; ++i) {
    if (c[i] == true) {
      lnl += log(pr0[i] + (pr1[i] * st[i])) - log(pr0[i] + (pr1[i] * st0[i]));
    } else {
      lnl += log(pr1[i]) + ln_ft[i]         - log(pr0[i] + (pr1[i] * st0[i]));
    }
  }
  
  return -lnl;
}')

cens <- (lo==0) | (rc==0)
fC <- spweib_lnlC(cens, pr0, pr1, st, st0, ln.ft)
fR <- cens5(cens, pr0, pr1, st, st0, ln.ft)

microbenchmark(
  cens5(cens, pr0, pr1, st, st0, ln.ft), 
  spweib_lnlC(cens, pr0, pr1, st, st0, ln.ft)
  )

sourceCpp("../core-functions/spweibull.cpp")


source("../core-functions/spweib-lnl-versions.r")


library(ggplot2)

foo <- microbenchmark(
  spweib_lnl(theta, Y, X, Z),
  spweib_lnl2(theta, Y, X, Z))

foo
autoplot(foo)

