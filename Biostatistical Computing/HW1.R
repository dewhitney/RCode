#Regression based on homoskedastic normal data of equal group sizes.
do.one = function(b0,b1,n){
  times = n/10
  x = sort(rep(seq(from=0.1,to=1.0,by=0.1),times))
  sims = function(x){return(rnorm(1,mean={b0+b1*x},sd=1))}
  y = mapply(x, FUN=sims)
  df = data.frame(x,y)
  model = with(df, lm(y~x))
  return(summary(model))
}

beta0 = c(0L,20L)
beta1 = c(0L,10L)
n = c(10L,100L,1000L)

#Had some problems working with the 3D array, so used same seed and stored
#estimated betas and std errors separately. Revisit this later. -DW
write(c("b0","b1","n","b0hat","b1hat","coverage"), file="HW1out.txt", append=FALSE, sep=",")
for (i in beta0){
  for(j in beta1){
    for (k in n){
      set.seed(1)
      many.betas = t(replicate(1000, do.one(b0=i,b1=j,n=k)$coefficients[,1]))
      many.mean = colSums(many.betas/1000)
      
      set.seed(1)
      many.SEs = t(replicate(1000, do.one(b0=i,b1=j,n=k)$coefficients[,2]))
      
      crit = qt(.025,df=k-1,lower.tail=FALSE)
      cover.beta1 = mean({{many.betas[,2] - crit*many.SEs[,2]} < j } & {j < {many.betas[,2] + crit*many.SEs[,2]}})
      write(c(i,j,k,as.vector(many.mean),cover.beta1), file="HW1out.txt", append=TRUE, sep=",")
    }
  }
}

#set.seed(1)
#many.betas = t(replicate(1000, do.one(beta0[1],beta1[1],n[1])$coefficients[,1]))
#many.mean = colSums(many.betas/1000)

#set.seed(1)
#many.SEs = t(replicate(1000, do.one(beta0[1],beta1[1],n[1])$coefficients[,2]))

#cover.beta1 = mean({{many.betas[,2] - many.SEs[,2]} < beta0[1] } & {beta0[1] < {many.betas[,2] + many.SEs[,2]}})

#Regression based on data with uniform errors of equal group sizes.
do.one = function(b0,b1,n){
  times = n/10
  x = sort(rep(seq(from=0.1,to=1.0,by=0.1),times))
  sims = function(x){return(b0+b1*x+runif(1,-3,3))}
  y = mapply(x, FUN=sims)
  df = data.frame(x,y)
  model = with(df, lm(y~x))
  return(summary(model))
}


write(c("b0","b1","n","b0hat","b1hat","coverage"), file="HW1out.txt", append=TRUE, sep=",")
for (i in beta0){
  for(j in beta1){
    for (k in n){
      set.seed(1)
      many.betas = t(replicate(1000, do.one(b0=i,b1=j,n=k)$coefficients[,1]))
      many.mean = colSums(many.betas/1000)
      
      set.seed(1)
      many.SEs = t(replicate(1000, do.one(b0=i,b1=j,n=k)$coefficients[,2]))
      
      crit = qt(.025,df=k-1,lower.tail=FALSE)
      cover.beta1 = mean({{many.betas[,2] - crit*many.SEs[,2]} < j } & {j < {many.betas[,2] + crit*many.SEs[,2]}})
      write(c(i,j,k,as.vector(many.mean),cover.beta1), file="HW1out.txt", append=TRUE, sep=",")
    }
  }
}

#Regression based on heteroskedastic normal data of equal group sizes.
do.one = function(b0,b1,n){
  times = n/10
  x = sort(rep(seq(from=0.1,to=1.0,by=0.1),times))
  sims = function(x){return(rnorm(1,mean={b0+b1*x},sd=x))}
  y = mapply(x, FUN=sims)
  df = data.frame(x,y)
  model = with(df, lm(y~x))
  return(summary(model))
}


write(c("b0","b1","n","b0hat","b1hat","coverage"), file="HW1out.txt", append=TRUE, sep=",")
for (i in beta0){
  for(j in beta1){
    for (k in n){
      set.seed(1)
      many.betas = t(replicate(1000, do.one(b0=i,b1=j,n=k)$coefficients[,1]))
      many.mean = colSums(many.betas/1000)
      
      set.seed(1)
      many.SEs = t(replicate(1000, do.one(b0=i,b1=j,n=k)$coefficients[,2]))
      
      crit = qt(.025,df=k-1,lower.tail=FALSE)
      cover.beta1 = mean({{many.betas[,2] - crit*many.SEs[,2]} < j } & {j < {many.betas[,2] + crit*many.SEs[,2]}})
      write(c(i,j,k,as.vector(many.mean),cover.beta1), file="HW1out.txt", append=TRUE, sep=",")
    }
  }
}
