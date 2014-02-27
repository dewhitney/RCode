## A. Implement the LASSO for linear and logistic regression.

# the gradient function for linear regression
linear = function(x,y,Beta){t(x)%*%(x%*%Beta-y)}

# the gradient function for logistic regression
expit = function(t){1/(1+exp(-t))}
p = function(x,Beta){expit(apply(x,1,function(u){t(u)%*%Beta}))}
logistic = function(x,y,Beta){-t(x)%*%(y-p(x,Beta))} #grad

# Requires soft thresholding
sft = function(w,lambda){
  apply(as.matrix(w),1,
        function(x){
          sign(x)*max((abs(x)-lambda),0)
          }
        )
}

## Lasso regression by gradient descent
# Inputs:
# x - a numeric data.frame, vector or matrix of predictor variables
# y - a numeric vector the outcome of interest
# grad - the gradient for the desired model (linear, logistic, etc)
# lambda - the penalty tuning parameter. Set=0 gives original glm estimates
# const - if TRUE, fits an intercept
# tol - stopping criteria for convergence
# Outputs:
# coef1,2,... - estimates for the coefficients in the model
# iterations - how many iterations occurred before estimates within tol
lasso.gd = function(x,y,grad,lambda,const=FALSE,tol=1e-3){
  x = as.matrix(x)
  if(const){x = cbind(1,x)} #determines presence of intercept
  k = ncol(x) #number of parameters
  Beta = numeric(k)
  t = 1/(max(svd(t(x)%*%x)$d))
  i = 0
  dBeta = 1
  while(dBeta > tol){
    Beta0 = Beta
    Beta = sft(Beta0 - t*grad(x,y,Beta0), t*lambda)
    dBeta = norm(Beta0 - Beta,"1")
    i = i+1
  }
  return(c(coef=Beta,iterations=i))
}

## Sanity check. Does it work?
set.seed(101990)
x = seq(0,1,1/9)
y = x+rnorm(10)
lm(y~0+x);lasso.gd(x,y,linear,lambda=0)

## Simulate normal data p>>n
set.seed(17)
n = 100
p = 200
results = data.frame(s=1:5,lambda=numeric(5))
not0betas = matrix(0,5,5)
lambda = 1
for(s in 1:5){
  x = matrix(rnorm(n*p),nrow=n,ncol=p)
  y = rowSums(x[,1:5]) + rnorm(n=n,sd=s)
  Betas = lasso.gd(x,y,linear,lambda,tol=1e-03)[-(p+1)]
  while(sum(abs(Betas[-201])>0)>5){
    lambda = lambda+1
    Betas = lasso.gd(x,y,linear,lambda,tol=1e-03)[-(p+1)]
  }
  results$lambda[s]=lambda
  not0betas[s,]=which(abs(Betas)>0)
}


library(reshape2)
## Randomly splitting my data into 10 groups
group=function(x){
  randomize = sample(0:99) %/% 10
  labeled = data.frame(x,randomize)
}

