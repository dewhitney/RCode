# David Whitney
# HW 3
# BIOST 562a
# Winter 2014
# Instructor: Noah Simon

# Problem 1

# First, a function for testing membership an element in an interval
is.in = function(x,L,U){L <= x & x <= U}

# This function calculates the p-value of the F-test for a model (lm object)
pval = function(model){
  pf(summary(model)$fstatistic[1],
     summary(model)$fstatistic[2],
     summary(model)$fstatistic[3],
     lower.tail=FALSE)
}

# This function creates a sequence of n evenly spaced observations x
# on the unit interval, then samples one y ~ N(mean = x, sd = 1)
# for each observation.
# Three models are fit:
# a) y ~ x (linear regression)
# b) y ~ 1st + 2nd + 3rd + 4th quartile (as categories/factors)
# c) y ~ 1st + 2nd + 3rd + 4th quartile (as numbers) (trend analysis)
# We return the result of the F-test for an association between
# y and the predictors in each case based on a 0.05 significance level.
# A result of TRUE means that we made a Type II error for that model.

do.typeII = function(n, probs = seq(0,1,0.25)){
  x = seq(0,1,by=1/(n-1))
  y1 = x + rnorm(n)
  y2 = 2.5*x^2 - 2*x + rnorm(n) 
  
  q = quantile(x, probs=probs)
  xbyq = character(n)
  for(i in 1:4){xbyq[is.in(x,q[i],q[i+1])]=i}
  
  lm1.a = lm(y1~x)
  lm1.b = lm(y1~xbyq)
  lm1.c = lm(y1~as.numeric(xbyq))
  
  lm2.a = lm(y2~x)
  lm2.b = lm(y2~xbyq)
  lm2.c = lm(y2~as.numeric(xbyq))
  
  alpha = 0.05
  p = c(A1=pval(lm1.a),B1=pval(lm1.b),C1=pval(lm1.c),
        A2=pval(lm2.a),B2=pval(lm2.b),C2=pval(lm2.c))
  return(p>alpha)
}

set.seed(19)
B = 10000
power = 1 - apply(replicate(B, do.typeII(n=100)),MARGIN=1,FUN=mean)
morefactor = 1 - apply(replicate(B, do.typeII(n=100,probs = seq(0,1,0.1))),MARGIN=1,FUN=mean)
lessfactor = 1 - apply(replicate(B, do.typeII(n=100,probs = seq(0,1,0.5))),MARGIN=1,FUN=mean)
allfactor = 1 - apply(replicate(B, do.typeII(n=100,probs = seq(0,1,0.01))),MARGIN=1,FUN=mean)

# Problem 2
# Coordinate descent

library(MASS)

lm.descent = function(x,y, epsilon=1e-10, const=TRUE){
  ols = coef(lm(y~x))
  if(const){x = cbind(1,x)}
  x = as.matrix(x)
  k = ncol(x)
  n = length(y)
  betas = numeric(k)
  not.done = TRUE
  ticker = 0
  while(not.done){
    #prev = betas
    for(i in 1:k){
      r = numeric(n)
      for(j in 1:n){
        r[j] = y[j] - t(betas[-k]) %*% x[j,-k]
      }
      betas[k] = sum(x[,k]*r) / (t(x[,k]) %*% x[,k])
    }
    not.done = t(ols-betas) %*% {ols-betas} > epsilon
    ticker = ticker+1
  }
  return( c(betas,after=ticker) )
}

set.seed(1990)
n = 100

xi = matrix(rnorm(2*n),ncol=2)
yi = x[,1] + rnorm(n)
lm.descent(xi,yi)

xii = mvrnorm(n, mu=c(0,0), Sigma = matrix(c(1,.8,.8,1),ncol=2))
yii = x[,1] + rnorm(n)
lm.descent(xii,yii)