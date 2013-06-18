# the msm library provides functions for truncated normal distributions
# install.packages("msm")
library(msm)

logL <- function(theta, x){
  return(sum(log(dtnorm(lower=theta[1], upper=theta[2], mean=theta[3], sd=theta[4], x[x>theta[1] & x<theta[2]]))))
}

logLxMinus <- function(theta, x=data){
  print(theta);
  erg <- (-1*logL(theta, x))
  print(erg)
  return(erg)
}

start <- c(1500, 17500, 9500, 1250)
lowerConst <- c(1000, 15000, 2000, 500)
upperConst <- c(2000, 20000, 15000, 2000)

data=c(rtnorm(100000,lower=1,mean=1,sd=1), rtnorm(100000,lower=1,mean=20,sd=3), rtnorm(10000,lower=1,mean=6000,sd=1000));

opt = optim(par=start, fn=logLxMinus, x=data, lower=lowerConst, upper=upperConst, method="L-BFGS-B")
plot(density(data[data>opt$par[1] & data<opt$par[2]]))
lines(seq(opt$par[1],opt$par[2]),dtnorm(lower=opt$par[1],upper=opt$par[2],mean=opt$par[3],sd=opt$par[4],seq(opt$par[1],opt$par[2])))

data_histo_raw = as.matrix(table(round(data)))
data_histo = cbind(as.numeric(rownames(data_histo_raw)), data_histo_raw[,1])
rownames(data_histo) = NULL

logLhisto <- function(theta, histo){
  histoInRange <- histo[histo[,1]>theta[1] & histo[,1]<theta[2],]
  return(sum(histoInRange[,2]*log(dtnorm(lower=theta[1], upper=theta[2], mean=theta[3], sd=theta[4], histoInRange[,1]))))
}

logLhistoMinus <- function(theta, histo){
  print(theta);
  erg <- (-1*logLhisto(theta, histo))
  print(erg)
  return(erg)
}

opt = optim(par=start, fn=logLhistoMinus, histo=data_histo, lower=lowerConst, upper=upperConst, method="L-BFGS-B")
data_histo_inrange = data_histo[data_histo[,1]>opt$par[1] & data_histo[,1]<opt$par[2],]
plot(x=data_histo_inrange[,1], y=data_histo_inrange[,2]/sum(data_histo_inrange[,2]))
lines(seq(opt$par[1],opt$par[2]),dtnorm(lower=opt$par[1],upper=opt$par[2],mean=opt$par[3],sd=opt$par[4],seq(opt$par[1],opt$par[2])))

plotResult <- function(result, histo){
  histo_inrange = histo[histo[,1]>result[1] & histo[,1]<result[2],]
  plot(x=histo_inrange[,1], y=histo_inrange[,2])
  lines(seq(result[1],result[2]),dtnorm(lower=result[1],upper=result[2],mean=result[3],sd=result[4],seq(result[1],result[2]))*sum(histo_inrange[,2]),col="blue",lwd=3)
}

plotResultFull <- function(result, histo){
  histo_inrange = histo[histo[,1]>result[1] & histo[,1]<result[2],]
  plot(x=histo[,1], y=histo[,2],ylim=c(0,max(histo_inrange[,2])))
  lines(seq(min(histo[,1]),max(histo[,1])),dtnorm(lower=result[1],upper=result[2],mean=result[3],sd=result[4],seq(min(histo[,1]),max(histo[,1])))*sum(histo_inrange[,2]),col="blue",lwd=3)
}

parameters=cbind(c(10000,15000,20000),c(55000,60000,65000),c(30000,40000,45000),c(3000,5000,10000))
rownames(parameters)=c("lowerConst","start","upperConst")
colnames(parameters)=c("lower","upper","mean","sd")

opt = optim(par=parameters[2,], fn=logLhistoMinus, histo=data_histo, lower=parameters[1,], upper=parameters[3,], method="L-BFGS-B")
plotResult(opt$par, data)