lm.coorDescent = function(x,y,tol=1e-10,const=TRUE){
  x = as.matrix(x)
  if(const){x = cbind(1,x)} #determines presence of intercept
  k = ncol(x) #number of parameters
  Beta = numeric(k)
  i = 0
  dBeta = 1
  while(dBeta > tol){
    Beta0 = Beta
    for(j in 1:k){
      r = y - x[,-j]%*%Beta[-j]
      Beta[j] = t(x[,j])%*%r/(t(x[,j])%*%x[,j])
    }
    dBeta = norm(Beta0 - Beta,"2")
    i = i+1
  }
  return(c(coef=Beta,iterations=i))
}

## Sample Usage
set.seed(1990)
n = 100

#First for some independent data
xi = matrix(rnorm(2*n),ncol=2)
yi = xi[,1] + rnorm(n)
lm.coorDescent(xi,yi)
coef(lm(yi~xi)) #Always good to verify results

#More iterations are required for correlated covariates
library(MASS)
xii = mvrnorm(n, mu=c(0,0), Sigma = matrix(c(1,.8,.8,1),ncol=2))
yii = xii[,1] + rnorm(n)
lm.coorDescent(xii,yii)
coef(lm(yii~xii)) #Always good to verify results

