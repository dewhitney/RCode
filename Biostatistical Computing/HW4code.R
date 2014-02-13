## David Whitney's Homework 4
## Biostatistics 562

### Problem 1

set.seed(1990)

person = function(genes=1000,de=100,mean){
  c(rnorm(genes-de),mean+rnorm(de))
}

### 1a.

experiment = function(N,mean1=0,mean2,genes=1000,de=100){
  class1 = replicate(N, person(mean=mean1))
  class2 = replicate(N, person(mean=mean2))
  return(cbind(class1,class2))
}

N = 50
XP = experiment(N=N,mean2=0.5)

class = 1:N
p.values = apply(XP,MARGIN=1,FUN=function(X){t.test(X[class],X[-class])$p.val})
adjust.vals = p.adjust(p.values,method="BH")

fdr = 0.10
n.discoveries = sum(adjust.vals < fdr)
n.false = sum(adjust.vals[1:900] < fdr)

### 1b.

do.one = function(N=50,meandiff=0.5,fdr=0.10){
  XP = experiment(N=N,mean2=meandiff)
  class = 1:N
  p.values = apply(XP,MARGIN=1,FUN=function(X){t.test(X[class],X[-class])$p.val})
  BH.vals = p.adjust(p.values,method="BH")
  bon.vals = p.adjust(p.values,method="bonferroni")
  
  BH.n.discoveries = sum(BH.vals < fdr)
  BH.n.true = sum(BH.vals[901:1000] < fdr)
  BH.n.false = sum(BH.vals[1:900] < fdr)
  BH.fdp = BH.n.false/BH.n.discoveries
  if(BH.n.discoveries==0){BH.fdp=0}
  BH.tdr = BH.n.true/100
  
  bon.n.discoveries = sum(bon.vals < fdr)
  bon.n.true = sum(bon.vals[901:1000] < fdr)
  bon.n.false = sum(bon.vals[1:900] < fdr)
  bon.fdp = bon.n.false/bon.n.discoveries
  if(bon.n.discoveries==0){bon.fdp=0}
  bon.tdr = bon.n.true/100  
  
  return(c(BH_FDP=BH.fdp,BH_TDR=BH.tdr,Bonf_FDP=bon.fdp,Bonf_TDR=bon.tdr))
}

set.seed(1776)
B = 1000

TrueFalse.1 = apply(replicate(B, do.one()), MARGIN=1, FUN=mean)

### Problem 2

### 2a.

set.seed(4951)
TrueFalse.2a = apply(replicate(B, do.one(meandiff=2)), MARGIN=1, FUN=mean)

### 2b.

set.seed(7756)
TrueFalse.2b = apply(replicate(B, do.one(meandiff=.2)), MARGIN=1, FUN=mean)

set.seed(4771)
means = seq(.01,1,by=.01)
TrueFalse.2c = apply(replicate(B, do.one(meandiff=means)), MARGIN=1, FUN=mean)

### Problem 3

person.t = function(genes=1000){
  c(rt(n=genes,df=3))
}

experiment.t = function(N){
  class1 = replicate(N, person.t())
  class2 = replicate(N, person.t())
  return(cbind(class1,class2))
}

do.one.t = function(N=10,meandiff=0,fdr=0.10){
  XP = experiment.t(N=N)
  class = 1:N
  p.values = apply(XP,MARGIN=1,FUN=function(X){t.test(X[class],X[-class])$p.val})
  BH.vals = p.adjust(p.values,method="BH")
  bon.vals = p.adjust(p.values,method="bonferroni")
  
  BH.n.discoveries = sum(BH.vals < fdr)
  
  bon.n.discoveries = sum(bon.vals < fdr)
  
  return(c(BH=BH.n.discoveries,Bonf=bon.n.discoveries))
}

TrueFalse.3 = apply(replicate(B, do.one.t()), MARGIN=1, FUN=mean)


do.one.z = function(N=10,meandiff=0,fdr=0.10){
  XP = experiment.t(N=N)
  class = 1:N
  p.values = apply(XP,MARGIN=1,FUN=function(X){pnorm(mean(X[class])-mean(X[-class]), mean=0, sd=sqrt(N*3))})
  BH.vals = p.adjust(p.values,method="BH")
  bon.vals = p.adjust(p.values,method="bonferroni")
  
  BH.n.discoveries = sum(BH.vals < fdr)
  
  bon.n.discoveries = sum(bon.vals < fdr)
  
  return(c(BH=BH.n.discoveries,Bonf=bon.n.discoveries))
}

TrueFalse.4 = apply(replicate(B, do.one.t()), MARGIN=1, FUN=mean)
