findMaxLowess <- function(dataX, lowerBound=0, upperBound=max(dataX[,1]), f=.1){
	xMax=0;
	while(xMax <= lowerBound+5 && upperBound > lowerBound+30){
		indicesInRange <- dataX[,1]<=upperBound & dataX[,1]>=lowerBound;
		maxLowess<-lowess(dataX[,1][indicesInRange], dataX[,2][indicesInRange], f=f);
		xMax<-which.max(maxLowess[[2]])[1];
		upperBound = lowerBound + 0.9*(upperBound-lowerBound);
	}
	return(maxLowess)
}
