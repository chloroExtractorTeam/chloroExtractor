findMaxLowess <- function(data){
	xMax = 0
	upperBound = max(data[,1])
	while(xMax < 1 & upperBound > 100){
		maxLowess <- lowess(data[,1][data[,1]<upperBound], data[,2][data[,1]<upperBound]);
		xMax <- which.max(maxLowess[[2]])[1];
		upperBound <- 0.9*upperBound;
	}
	return(maxLowess)
}