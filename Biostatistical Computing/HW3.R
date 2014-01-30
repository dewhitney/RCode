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
  y = x + rnorm(n)
  
  q = quantile(x, probs=probs)
  xbyq = character(n)
  for(i in 1:4){xbyq[is.in(x,q[i],q[i+1])]=i}
  
  lm.a = lm(y~x)
  lm.b = lm(y~xbyq)
  lm.c = lm(y~as.numeric(xbyq))
  
  alpha = 0.05
  p = c(A=pval(lm.a),B=pval(lm.b),C=pval(lm.c))
  return(p<alpha)
}

set.seed(19)
B = 10000
power = apply(replicate(B, do.typeII(n=n)),MARGIN=1,FUN=mean)
