#' Regular Weibull regression
weibullCpp <- function(Y, X, inits=NULL, max.iter, silent=TRUE) {
  if (is.null(inits)) {
    inits <- c(rep(0, ncol(X)), 0)
  }
  
  trace <- !silent
  est <- optim(inits, weib_lnlCpp, method="BFGS", 
               control=list(trace=trace, maxit=max.iter), 
               hessian=T, y=Y, X=X)
  
  # Solve other results
  if (est$convergence!=0 & !silent) stop('Model did not converge')
  coef <- est$par
  vcv <- solve(est$hessian)
  vcv <- make.positive.definite(vcv)
  logL <- -est$value
  
  # Put together results
  return(list(coefficients = coef, vcv = vcv, logL = logL))
}

# Split Population
spweibullCpp <- function(Y, X, Z, max.iter, silent=FALSE) {  
  # Estimate base model
  if (!exists("base.inits")) {
    base.inits <- c(rep(0, ncol(X)), 0)
  }
  if (!silent) cat('Fitting base weibull...\n')
  base <- weibull(Y=Y, X=X, inits=base.inits, max.iter=200, silent=TRUE)
  
  # Estimate full model
  x.inits <- base$coefficients[1:ncol(X)]
  a.init <- base$coefficients[ncol(X)+1]
  if (!silent) cat('Fitting split weibull...\n')
  trace <- !silent
  est <- optim(c(x.inits, rep(0, ncol(Z)), a.init), spweib_lnlCpp, method="BFGS", 
    control=list(trace=trace, maxit=max.iter), hessian=T, y=Y, X=X, Z=Z)
  
  # Solve other results
  if (est$convergence!=0 & !silent) stop('Model did not converge')
  coef <- est$par
  vcv <- solve(est$hessian)
  vcv <- make.positive.definite(vcv)
  logL <- -est$value
  
  # Put together results
  return(list(coefficients = coef, vcv = vcv, logL = logL, base=base))
}

#' Regular Log-logistic regression
loglogCpp <- function(Y, X, inits=NULL, max.iter, silent=TRUE) {
  if (is.null(inits)) {
    inits <- c(rep(0, ncol(X)), 0)
  }
  
  trace <- !silent
  est <- optim(inits, loglog_lnlCpp, method="BFGS", 
               control=list(trace=trace, maxit=max.iter), 
               hessian=T, y=Y, X=X)
  
  # Solve other results
  if (est$convergence!=0 & !silent) stop('Model did not converge')
  coef <- est$par
  vcv <- solve(est$hessian)
  vcv <- make.positive.definite(vcv)
  logL <- -est$value
  
  # Put together results
  return(list(coefficients = coef, vcv = vcv, logL = logL))
}

# Split Population
sploglogCpp <- function(Y, X, Z, max.iter, silent=FALSE) {
  # Estimate base model
  if (!exists("base.inits")) {
    base.inits <- c(rep(0, ncol(X)), 1)
  }
  if (!silent) cat('Fitting base loglog...\n')
  base <- loglog(Y=Y, X=X, inits=base.inits, max.iter=200, silent=TRUE)
    
  # Estimate full model
  x.inits <- base$coefficients[1:ncol(X)]
  a.init <- base$coefficients[ncol(X)+1]
  if (!silent) cat('Fitting split loglog...\n')
  trace <- !silent
  est <- optim(c(x.inits, rep(0, ncol(Z)), a.init), sploglog_lnlCpp, method="BFGS", 
    control=list(trace=trace, maxit=max.iter), hessian=T, y=Y, X=X, Z=Z)
  
  # Solve other results
  if (est$convergence!=0 & !silent) stop('Model did not converge')
  coef <- est$par
  vcv <- solve(est$hessian)
  vcv <- make.positive.definite(vcv)
  logL <- -est$value
  
  # Put together results
  return(list(coefficients = coef, vcv = vcv, logL = logL, base=base))
}

spdurCpp <- function(duration, atrisk, data=NULL, last="end.spell", t.0="t.0", 
                  fail="failure", distr='weibull', max.iter=300, na.action, 
                  silent=FALSE, ...) 
{ 
  cl <- match.call()
  
  # NA actions, separately because model.frame for two different equations
  # might return data frames with different row numbers, breaking the function.
  f1 <- as.formula(eval(duration))
  f2 <- as.formula(eval(atrisk))
  vars <- unique(c(all.vars(f1), all.vars(f2)))
  vars <- c(vars, last, t.0, fail)
  if (missing(na.action)) na.action <- options("na.action")[[1]]
  df <- do.call(na.action, list(data[, vars]))
  
  # Duration equation
  mf.dur <- eval(model.frame(formula=duration, data=df), parent.frame())
  X <- model.matrix(attr(mf.dur, 'terms'), data=mf.dur)
  attr(X, "na.action") <- na.action(df)
  attr(X, "dimnames")[[2]] <- sub("\\(Intercept\\)", "(Dur. Intercept)", 
                                  attr(X, "dimnames")[[2]])
  lhb <- model.response(mf.dur)
  # Risk/non-immunity equation
  mf.risk <- eval(model.frame(formula=atrisk, data=df), parent.frame())
  Z <- model.matrix(attr(mf.risk, 'terms'), data=mf.risk)
  attr(Z, "na.action") <- na.action(df)
  attr(Z, "dimnames")[[2]] <- sub("\\(Intercept\\)", "(Risk Intercept)", 
                                  attr(Z, "dimnames")[[2]])
  lhg <- model.response(mf.risk)
  # Y vectors
  Y <- cbind(atrisk=lhg, duration=lhb, last=df[, last], t.0=df[, t.0], 
             fail=df[, fail])
  attr(Y, "last") <- last
  attr(Y, "t.0") <- t.0
  attr(Y, "fail") <- fail
  
  # Estimation
  if (distr=='weibull') {
    est <- spweibullCpp(Y, X, Z, max.iter, silent=silent, ...)
  }
  if (distr=='loglog') {
    est <- sploglogCpp(Y, X, Z, max.iter, silent=silent, ...)
  }
  # Names
  varnames <- c(paste(unlist(attr(X, 'dimnames')[2])), paste(unlist(attr(Z, 'dimnames')[2])), 'log(alpha)')
  attr(est$coefficients, 'names') <- varnames
  colnames(est$vcv) <- rownames(est$vcv) <- varnames
  
  # Calculate uncertainty measures
  est$se <- sqrt(diag(est$vcv))
  est$zstat <- with(est, coefficients/se)
  est$pval <- 2*(1-pnorm(abs(est$zstat)))
  
  # Other class elements
  est$mf.dur <- mf.dur
  est$mf.risk <- mf.risk
  est$Y <- Y
  est$na.action <- attr(df, "na.action")
  est$call <- cl
  est$distr <- eval.parent(distr, 1)
  est$obs <- nrow(Y)
  
  class(est) <- 'spdur'
  return(est)
}
