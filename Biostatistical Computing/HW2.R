# Part 1
# Determine how many significant figures to report coverage
# probabilities to.
qbinom(p=c(.05,.95), size=1000, prob=.95)


# Part 2
library(sandwich) #needed for vcovHC() robust std errors.
library(MASS) #needed for ginv() matrix inverses.

# This function takes an lm object and calculates robust confidence intervals
# it is based on the guts of the confint.default function
ci.robust = function(object,parm,level = 0.95,...)
{
  cf = coef(object)
  pnames = names(cf)
  if (missing(parm))
    parm = pnames
  else if (is.numeric(parm))
    parm = pnames[parm]
  a = (1 - level)/2
  a = c(a,1-a)
  pct = round(100*a,1)
  fac = qnorm(a)
  ci = cf + sqrt(diag( vcovHC(object, "HC0") )) %o% qnorm(a)
  ci
}

# Part 3
# Code from Noah Simon's HW Key 1:
# a little helper function to save typing later
is.element <- function(x, interval){
  (x > min(interval)) & (x < max(interval))
}

do.one <- function(beta0, beta1, n){
  x <- rep(1:10 /10, n/10)
  mean.y <- beta0 + beta1*x
  y3 <- mean.y + rnorm(n, mean=0, sd=sqrt(x))
  lm3 <- lm(y3~x)
  r.ci3 <- ci.robust(lm3)
  ci3 <- confint.default(lm3)
  return(c(
    betahat.3 = coef(lm3), cov.3 = is.element(beta1, ci3[2,]), rocov.3 = is.element(beta1, r.ci3[2,])
  ))
}

# Note that we can get expected values of betahat AND the coverage
# by taking the mean() of rows of replicates
set.seed(1990)
B = 10000
a1 <- apply( replicate(B, do.one(0,1, 10) ) , 1, mean)
a2 <- apply( replicate(B, do.one(0,1, 100) ) , 1, mean)
a3 <- apply( replicate(B, do.one(0,1,1000) ) , 1, mean)

# Part 4
# The following function does matrix multiplications
# to obtain point estimates for a simple regression and then
# estimates their robust standard errors.
simple.robust.lm = function(x,y){
  n = length(x)
  X = matrix( c(rep(1,n),x), nrow=n, ncol=2)
  H = ginv(t(X) %*% X) %*% t(X)
  betas = H %*% y
  error = y - X %*% betas
  cov = H %*% diag(c(error^2)) %*% t(H)
  return(c(b0=betas[1],b1=betas[2],cov0=diag(cov)[1],cov1=diag(cov)[2]))
}

# This function calculates confidence intervals
# but cannot take class = lm arguments.
# model[1:2] must contain coefficient estimates
# model[3:4] must contain (robust) std error estimates.
ci.robust.matrix = function(model, level=.95)
{
  a = (1 - level)/2
  a = c(a,1-a)
  ci = model[1:2] + sqrt(model[3:4]) %o% qnorm(a)
  ci
}
