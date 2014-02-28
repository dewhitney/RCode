# Soft thresholding for penalty step in gradient descent
sft = function(w,lambda){
	apply(as.matrix(w),1,
		function(x){sign(x)*max((abs(x)-lambda),0)}
		)
}

# the gradient function for linear regression
linear = function(x,y,Beta){t(x)%*%(x%*%Beta-y)}

# the gradient function for logistic regression
expit = function(t){1/(1+exp(-t))}
p = function(x,Beta){expit(apply(x,1,function(u){t(u)%*%Beta}))}
logistic = function(x,y,Beta){-t(x)%*%(y-p(x,Beta))} #grad

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
# nonZero - how many non-zero coefficients are in the final model
lasso.gd = function(x,y,grad=linear,lambda,const=FALSE,tol=1e-3){
	x = as.matrix(x)
	if(const){x = cbind(1,x)} #determines presence of intercept
	k = ncol(x) #number of parameters
	Beta = numeric(k); t = 1/(max(svd(t(x)%*%x)$d)); i = 0; dBeta = 1
	while(dBeta > tol){
		Beta0 = Beta
		if(const==TRUE){
			Beta[1]=Beta0[1] - t*grad(x,y,Beta0)[1]
			Beta[-1]=sft(Beta0[-1] - t*grad(x,y,Beta0)[-1], t*lambda)
			}
		else{
			Beta = sft(Beta0 - t*grad(x,y,Beta0), t*lambda)
			}
		dBeta = norm(Beta0 - Beta,"2")
		i = i+1
	}
	nBetas = sum(Beta>0)
	return( list(coefficients=Beta,extra=c(iterations=i,nonZero=nBetas)) )
}

# Simulate data for experiments
rData = function(sigma=1,n=100,p=200){
	x = matrix(rnorm(n*p),nrow=n,ncol=p)
	y = rowSums(x[,1:5]) + rnorm(n=n,sd=sigma)
	return(data.frame(y,x))
}

# Randomly splitting my data into 10 groups
tenfolds=function(x){
	n = dim(x)[1]-1
	grp = sample(0:n) %/% ((n+1)/10)
	labeled = data.frame(grp,x)
	return(labeled)
}

# Perform cross validation on ten folds of the data
cv10f.lasso=function(x,lambda){
	obs=tenfolds(x)
	SS = numeric(10)
	for(i in 0:9){
		train=lasso.gd(obs[which(obs$grp != i),3:202],obs[which(obs$grp != i),2],linear,lambda)$coefficients
		fit=as.matrix(obs[which(obs$grp == i),3:202]) %*% as.matrix(train)
		y=as.matrix(obs[which(obs$grp == i),2])
		SS[i+1]=sum((y-fit)^2)
		}
	return(c(Lambda=lambda,SSR=sum(SS)))
}

# Find the best lambda on a grid
best.lambda=function(Data,lower,upper,delta){
	grid = seq(lower,upper,delta)
	m = length(grid)
	results=data.frame(Lambda=grid,SSR=numeric(m))
	i = 0
	for(lambda in grid){
		i = i+1 
		results[i,] = cv10f.lasso(Data,lambda)
	}
	bestLambda=with(results,Lambda[which.min(SSR)])
	bestFit=lasso.gd(Data[,-1],Data[,1],grad=linear,lambda=bestLambda)
	return(c(lambda=bestLambda,non.zero=bestFit$extra[2]))
}

do.one=function(s=1,l,u,d){
	Data=rData(sigma=s)
	return(best.lambda(Data,l,u,d))
	}

#######################################################################
#######################################################################

set.seed(509)

answer = data.frame(s=1:5,lambda=0,correct=0)

for(s in 1:5){
	vars = 6
	Data=rData(sigma=s)
	lambda = 1
	while(vars > 5){
		lambda = lambda+1
		fit = lasso.gd(Data[,-1],Data[,1],grad=linear,lambda)$coefficients
		vars = sum(fit>0)
	}
	correct = sum(fit[1:5]>0)
	answer$lambda[s] = lambda
	answer$correct[s] = correct
}

do.one(s=1,l=3,u=21,d=3)

B=10
meanResults = data.frame(Lambda=numeric(3), Variables=numeric(3)) 
for(i in 1:3){
	meanResults[i,]=apply(replicate(B,do.one(s=i,l=3,u=99,d=3)),1,mean)
	}

