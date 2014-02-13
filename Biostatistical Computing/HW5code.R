## David Whitney
## Biostatistical Computing
## Winter 2014
## Homework 5

#########################################################################
#########################################################################

## Problem 1

# First, a function for testing membership an element in an interval
is.in = function(x,L,U){L <= x & x <= U}

# Two tests for the multiple chisquare and f tests of glm or lm objects
glm.test = function(model){anova(model,test="Chisq")["Pr(>Chi)"][2,1]}
lm.test = function(model){
  pf(
    q=summary(model)$fstatistic[1],
    df1=summary(model)$fstatistic[2],
    df2=summary(model)$fstatistic[3],
    lower.tail=FALSE)
}

# This function is for generating our random sample
data = function(n){
  x = seq(from=0,to=1,by=1/(n-1))
  q = quantile(x)
  xbyq = character(n)
  for(i in 1:4){xbyq[is.in(x,q[i],q[i+1])]=i}
  
  p = 1/(1+exp(-x))
  y = rbinom(n,prob=p,size=1) #y_i ~ bernoulli(p_i) with p_i=logit(x_i)
  z = rbinom(n,prob=.5,size=1)
  return(data.frame(y,z,x,xbyq))
}

do.one = function(n=100){
  experiment = data(n)
  logisticA = glm(y~x,family=binomial,data=experiment)
  logisticB = glm(y~xbyq,family=binomial,data=experiment)
  linearC = lm(y~x,data=experiment)
  linearD = lm(y~xbyq,data=experiment)
  
  ## using a level 0.05 test, compare the Type II error of approaches
  pvals1 = c(A=glm.test(logisticA),
            B=glm.test(logisticB),
            C=lm.test(linearC),
            D=lm.test(linearD))
  TypeII = pvals1 > 0.05
  
  
  logisticA = glm(z~x,family=binomial,data=experiment)
  logisticB = glm(z~xbyq,family=binomial,data=experiment)
  linearC = lm(z~x,data=experiment)
  linearD = lm(z~xbyq,data=experiment)
  
  ## using a level 0.05 test, compare the Type I error of approaches
  pvals2 = c(A=glm.test(logisticA),
             B=glm.test(logisticB),
             C=lm.test(linearC),
             D=lm.test(linearD))
  TypeI = pvals2 < 0.05
  
  return(c(TypeII=TypeII,TypeI=TypeI))
}

set.seed(509)
B = 10000
#apply(replicate(B,do.one()),1,mean) # Results shown below. Rerun takes TIME

## TypeII.A   TypeII.B   TypeII.C   TypeII.D        
## 0.7058     0.8152     0.7129     0.8264          
## TypeI.A    TypeI.B    TypeI.C    TypeI.D 
## 0.0522     0.0593     0.0503     0.0581 

#########################################################################
#########################################################################

## Problem 2

expit = function(t){1/(1+exp(-t))}
p = function(x,Beta){expit(apply(x,1,function(u){t(u)%*%Beta}))}
grad = function(x,y,Beta){-t(x)%*%(y-p(x,Beta))}

# Gradient Descent Implementation
# If I change the grad function above, this is a general optimizer
lm.gradescent = function(x,y,tol=1e-10,t=.05,const=TRUE){
  x = as.matrix(x)
  if(const){x = cbind(1,x)} #determines presence of intercept
  k = ncol(x) #number of parameters
  Beta = numeric(k)
  i = 0
  dBeta = 1
  while(dBeta > tol){
    Beta0 = Beta
    Beta = Beta0 - t*grad(x,y,Beta0)
    dBeta = norm(Beta0 - Beta,"2")
    i = i+1
  }
  return(c(coef=Beta,iterations=i))
}

## Does it work? (Yes.)
test = data(100)
coef(with(test, glm(y~x,family=binomial)))
with(test, lm.gradescent(x,y,t=.05))

# Newton-Raphson Implementation

# Requires second partial derivatives
grad2 = function(x,y,Beta){
  t(x)%*% diag(p(x,Beta)*(1-p(x,Beta))) %*%x
}

lm.NR = function(x,y,tol=1e-10,const=TRUE){
  x = as.matrix(x)
  if(const){x = cbind(1,x)} #determines presence of intercept
  k = ncol(x) #number of parameters
  Beta = numeric(k)
  i = 0
  dBeta = 1
  while(dBeta > tol){
    Beta0 = Beta
    Beta = Beta0 - ginv(grad2(x,y,Beta0))%*%grad(x,y,Beta0)
    dBeta = norm(Beta0 - Beta,"2")
    i = i+1
  }
  return(c(coef=Beta,iterations=i))
}

## Does it work? (Yes.)
test = data(100)
coef(with(test, glm(y~x,family=binomial)))
with(test, lm.NR(x,y))

## Comparisons of Gradient Descent and Newton-Raphson

set.seed(1990)
n = 100

x_a = matrix(rnorm(2*n),ncol=2)
y_a = rbinom(n,prob=expit(x_a[,1]),size=1)
lm.gradescent(x_a,y_a) # 28 iterations
lm.NR(x_a,y_a) # 6 iterations

x_b = mvrnorm(n, mu=c(0,0), Sigma = matrix(c(1,.8,.8,1),ncol=2))
y_b = rbinom(n,prob=expit(x_a[,1]),size=1)
lm.gradescent(x_b,y_b,t=.01) # 505 iterations, had to tune t
lm.NR(x_b,y_b) # 5 iterations

#########################################################################
#########################################################################

## Problem 3
## Try different numbers of covariates for independent and dependent data.

do.iterations = function(p,n=100,Cov=.8){
  x_a = matrix(rnorm(p*n),ncol=p)
  y_a = rbinom(n,prob=expit(x_a[,1]),size=1)
  gd_a = lm.gradescent(x_a,y_a,t=.01)[p+2]
  nr_a = lm.NR(x_a,y_a)[p+2]
  
  CovM = matrix(Cov,nrow=p,ncol=p)
  diag(CovM) = 1
  
  x_b = mvrnorm(n, mu=numeric(p), Sigma = CovM)
  y_b = rbinom(n,prob=expit(x_a[,1]),size=1)
  gd_b = lm.gradescent(x_b,y_b,t=.01)[p+2]
  nr_b = lm.NR(x_b,y_b)[p+2]
  return(c(GD.ind=gd_a,NR.ind=nr_a,GD.dep=gd_b,NR.dep=nr_b))
}

set.seed(1990)
P = 1:10
Results = apply(as.matrix(P),1,do.iterations)

##   Results
##                    [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
## GD.ind.iterations  131  136  192  159  205  194  213  189  318   456
## NR.ind.iterations    5    5    6    6    6    6    6    6    7     7
## GD.dep.iterations  117  317  575  500  504  533  557  642  612   680
## NR.dep.iterations    5    4    5    4    5    5    5    5    5     5
