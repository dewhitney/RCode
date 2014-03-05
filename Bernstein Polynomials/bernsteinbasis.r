## Bernstein polynomials
library(plyr)

# Mass function for truncated geometric distribution
dgeom.n=function(x, size, prob){
  d=dgeom(x, prob, log=FALSE)/pgeom(size, prob, lower.tail=TRUE, log.p=FALSE)
  return(d)
}

mkBnk = function(f=function(x){1},n,k,pmf=dbinom){
  coef = f(k/n)
  Bnk = function(x){coef*pmf(k,n,x)}
  return(Bnk)
}

# Expansion for a function
Bn = function(f=function(x){1},n,pmf=dbinom){
  g = function(x){f(x)}
  k = matrix(0:n)
  bnk = apply(k,1,FUN=function(m){mkBnk(f=g,k=m,n=n,pmf=pmf)})
  B = function(t){sum(sapply(t, each(bnk)))}
  return(function(t){apply(as.matrix(t),1,B)}) #the default f(x)=1 will satisfy B(X)=1
}

## de Casteljau's Algorithm for evaluating Bezier curves

deCasteljau = function(f,t,n){
  pts = seq(0,1,1/n)
  t0 = c(1-t,t)
  Beta = f(pts)
  for (k in 1:n){
    Beta = cbind(Beta[-(length(Beta))],Beta[-1]) %*% t0
    #print(Beta)
  }
  return(Beta)
}

fastB = function(f,n){
  B = function(t){return(deCasteljau(f=f,t,n=n))}
  return(function(t){apply(as.matrix(t),1,B)})
}

## N-th Forward difference function
## debug this
diff = function(f,x,dx,N=1,fwd=TRUE){
  pts=ifelse(fwd,seq(N*dx,0,-dx),seq(0,-N*dx,-dx))
  fpts=f(x+pts)
  terms=0:N
  coef=((-1)^terms)*choose(N,terms)
  return(cbind(fpts,coef))
}

#diff(function(x){x},0,.01)
