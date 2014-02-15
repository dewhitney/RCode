## Bernstein polynomials
library(plyr)

# Make the (n,k)th Bernstein Basis Polynomial
mkBnk = function(f=function(x){1},n,k){
  coef = choose(n,k)*f(k/n)
  p1 = k
  p2 = n-k
  Bnk = function(x){coef*x^p1*(1-x)^p2}
  #return(c(coef,p1,p2))
  return(Bnk)
}

# Expansion for a function
Bn = function(f=function(x){1},n){
  g = function(x){f(x)}
  k = matrix(0:n)
  bnk = apply(k,1,FUN=function(m){mkBnk(f=g,k=m,n=n)})
  B = function(t){sum(sapply(t, each(bnk)))}
  return(function(t){apply(as.matrix(t),1,B)}) #the default f(x)=1 will satisfy B(X)=1
}

## Example usage
#n = 2
#B20 = mkBnk(n=n,k=0); B21 = mkBnk(n=n,k=1); B22 = mkBnk(n=n,k=2)
#
#t = seq(0,1,.1)
#y20 = B20(t); y21 = B21(t); y22 = B22(t)
#plot(x=t,y=y20,type="l",lwd=2,col="blue")
#points(x=t,y=y21,type="l",lwd=2,col="forestgreen")
#points(x=t,y=y22,type="l",lwd=2,col="red")



## To use the output function, BX, on a vector of values t

#t = seq(0,1,.01)
#Sin = function(x){sin(2*pi*x)}
#BSin1 = Bn(f=Sin,n=1)
#BSin10 = Bn(f=Sin,n=10)
#BSin100 = Bn(f=Sin,n=100)
#fit1 = BSin1(t)
#fit10 = BSin10(t)
#fit100 = BSin100(t)
#plot(t,Sin(t),type="l",lwd=3,col="lightgreen",main="Approximating Sin by Bernstein Expansion")
#points(t,Sin(t),type="l",lwd=3,col="forestgreen")
#points(t,fit1,type="l",lwd=2,col="red")
#points(t,fit10,type="l",lwd=2,col="purple")
#points(t,fit100,type="l",lwd=2,col="blue")

