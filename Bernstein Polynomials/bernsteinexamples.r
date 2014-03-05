# Examples for Bernstein Polynomials

## Nice convergence for a concave function
f = function(x){-2*(x-.5)^2+1}
x = seq(0,1,.01)
plot(x,f(x),type="l")
for(N in 1:30){points(x,Bn(f,n=N,pmf=dbinom)(x),type="l")}

## Example with Sin scaled to interval 0 < x < 1.
t = seq(0,1,.01)
Sin = function(x){sin(2*pi*x)}
BSin1 = Bn(f=Sin,n=1)
BSin10 = Bn(f=Sin,n=10)
BSin100 = Bn(f=Sin,n=100)
fit1 = BSin1(t)
fit10 = BSin10(t)
fit100 = BSin100(t)
plot(t,Sin(t),type="l",lwd=3,col="lightgreen",main="Approximating Sin by Bernstein Expansion")
points(t,Sin(t),type="l",lwd=3,col="forestgreen")
points(t,fit1,type="l",lwd=2,col="red")
points(t,fit10,type="l",lwd=2,col="purple")
points(t,fit100,type="l",lwd=2,col="blue")

## Using de Casteljau's Algorithm
deCasteljau(Sin,.6,10)

fastBSin100 = fastB(Sin,100)

## Speed comparison
library(microbenchmark)
microbenchmark(BSin100(t),fastBSin100(t))
## Unit: milliseconds
## expr       min        lq    median       uq       max neval
## BSin100(t) 1213.9194 1233.5991 1247.0700 1289.011 1944.6590   100
## fastBSin100(t)  309.1887  316.9562  322.6283  342.625  564.1454   100