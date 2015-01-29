spweib_lnl <- function(theta, y, X, Z) {
  # Parameters
  rx <- ncol(X)
  rz <- ncol(Z)
  beta <- theta[1:rx]
  gamma <- theta[(rx+1):(rx+rz)]
  p <- theta[rx+rz+1]
  lambda <- exp(-X%*%beta)
  alpha <- exp(-p)
  
  # Dependent variables
  rc <- y[,1]    # right-censored spell (1 if uncensored)
  ti <- y[,2]
  t0 <- y[,4]    # duration at previous period
  lo <- y[,3]    # last observation in spell
  
  # Logistic for population split (1 = at risk)
  pr1 <- plogis(Z%*%gamma)
  pr0 <- plogis(Z%*%gamma, lower.tail=FALSE)
  
  # Weibull 
  #ln.ft <- log(alpha * lambda) + (alpha-1) * log(lambda * ti) - (lambda * ti)^alpha
  ln.ft <- log(alpha) + (alpha) * log(lambda) + (alpha - 1) * log(ti) - (lambda * ti)^alpha
  st    <- exp( -(lambda * ti)^alpha )    # S(t) for current period
  st0   <- exp( -(lambda * t0)^alpha )    # S(t) for previous period
  
  # Likelihood contributions
  cens   <- ifelse((lo==0) | (rc==0), log(pr0 + (pr1 * st)) - log(pr0 + (pr1 * st0)), 0)
  nocens <- ifelse((rc==1) & (lo==1), log(pr1) + ln.ft      - log(pr0 + (pr1 * st0)), 0)
  
  logl <- sum(cens+nocens)
  return(-logl)
}

spweib_lnl2 <- function(theta, y, X, Z) {
  # last ifelse() lnl calculation switched to Cpp
  # Parameters
  rx <- ncol(X)
  rz <- ncol(Z)
  beta <- theta[1:rx]
  gamma <- theta[(rx+1):(rx+rz)]
  p <- theta[rx+rz+1]
  lambda <- exp(-X%*%beta)
  alpha <- exp(-p)
  
  # Dependent variables
  rc <- y[,1]    # right-censored spell (1 if uncensored)
  ti <- y[,2]
  t0 <- y[,4]    # duration at previous period
  lo <- y[,3]    # last observation in spell
  
  # Logistic for population split (1 = at risk)
  pr1 <- plogis(Z%*%gamma)
  pr0 <- plogis(Z%*%gamma, lower.tail=FALSE)
  
  # Weibull 
  #ln.ft <- log(alpha * lambda) + (alpha-1) * log(lambda * ti) - (lambda * ti)^alpha
  ln.ft <- log(alpha) + (alpha) * log(lambda) + (alpha - 1) * log(ti) - (lambda * ti)^alpha
  st    <- exp( -(lambda * ti)^alpha )    # S(t) for current period
  st0   <- exp( -(lambda * t0)^alpha )    # S(t) for previous period
  
  # Likelihood contributions
  cens <- (lo==0) | (rc==0)
  logl <- spweib_lnlC(cens, pr0, pr1, st, st0, ln.ft)
  
  logl
}