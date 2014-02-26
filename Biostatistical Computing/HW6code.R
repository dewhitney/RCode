## A. Implement the LASSO for linear and logistic regression.

# the gradient function for linear regression
linear = function(x,y,Beta){-t(x)%*%(x%*%Beta-y)/length(y)}

# the gradient function for logistic regression
expit = function(t){1/(1+exp(-t))}
p = function(x,Beta){expit(apply(x,1,function(u){t(u)%*%Beta}))}
logistic = function(x,y,Beta){-t(x)%*%(y-p(x,Beta))} #grad

# Requires soft thresholding
sft = function(w,lambda){
  apply(as.matrix(w),1,function(x){sign(x)*(abs(x)-lambda)})
}

lasso.gradescent = function(x,y,grad,lambda,const=TRUE,tol=1e-10){
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
    dBeta = norm(Beta0 - Beta,"2")
    i = i+1
  }
  return(c(coef=Beta,iterations=i))
}

## Does it work?
with(test, lasso.gradescent(x,y,logistic,lambda=0,const=TRUE,tol=1e-10))


