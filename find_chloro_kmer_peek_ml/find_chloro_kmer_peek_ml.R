library(msm)

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
