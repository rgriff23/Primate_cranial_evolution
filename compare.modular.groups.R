
# Function to test for a significant difference between two groups of individuals in the strength of modularity 

groups = dimnames(A)[[3]]
groups[c(1,2,44:63)] = "Prosimian"
groups[3:43] = "Anthropoid"

groups2 = dimnames(A)[[3]]
groups2[1:length(groups2)] = "NonCat"
groups2[4:30] = "Cat"

compare.modular.groups <- function (A, modules, groups, ShowPlot=FALSE, iter=999) {
	
	# preparations
	groups <- as.factor(groups)
	modules <- as.factor(modules)
	nmod <- nlevels(modules)	
	if (nlevels(groups) != 2) {stop("Function only defined for two groups.")}
	if(any(is.na(A)) == T) {stop("Data matrix contains missing values.")}

	# define two arrays
	if (length(dim(A)) == 3) {
		if (length(modules) != dim(A)[1]) {stop("Not all landmarks are assigned to a partition.")}
		if (length(groups) != dim(A)[3]) {stop("Not all individuals are assigned to a group.")}
		A1 <- A[,,which(groups == levels(groups)[1])]
		A1 <- gpagen(A1, ShowPlot=FALSE)[[1]]
		A2 <- A[,,which(groups == levels(groups)[2])]
		A2 <- gpagen(A2, ShowPlot=FALSE)[[1]]
		x1 <- two.d.array(A1)
		x2 <- two.d.array(A2)
        p <- dim(A1)[1]
        k <- dim(A1)[2]
        gps <- as.factor(rep(modules, k, each = k, length=p*k))}
  	if (length(dim(A))==2){
        if(length(modules) != ncol(x)) {stop("Not all variables are assigned to a partition.")}
        if(length(groups) != nrow(x)) {stop("Not all individuals are assigned to a group.")}
        x1 <- A[which(groups == levels(groups)[1]),]
        x2 <- A[which(groups == levels(groups)[2]),]
        gps <- modules}
	
	# compute RV1
	S1 <- cov(x1)
  	RV1 <- array(0, dim=c(nmod, nmod))
 	for (i in 1:(nmod-1)){
    	for (j in 2: nmod){
     		S11 <- S1[which(gps == levels(gps)[i]), which(gps == levels(gps)[i])]
      		S22 <- S1[which(gps == levels(gps)[j]), which(gps == levels(gps)[j])]
      		S12 <- S1[which(gps == levels(gps)[i]), which(gps == levels(gps)[j])]
      		S21 <- t(S12)
      		RV1[i,j] <- sum(diag(S12%*%S21))/sqrt(sum(diag(S11%*%S11))*sum(diag(S22%*%S22)))
      		diag(RV1) <- 0
    	}
 	 }
 	 RV1 <- sum(RV1)/(nmod/2*(nmod-1))
	
	# compute RV2
	S2 <- cov(x2)
  	RV2 <- array(0, dim=c(nmod, nmod))
 	for (i in 1:(nmod-1)){
    	for (j in 2: nmod){
     		S11 <- S2[which(gps == levels(gps)[i]), which(gps == levels(gps)[i])]
      		S22 <- S2[which(gps == levels(gps)[j]), which(gps == levels(gps)[j])]
      		S12 <- S2[which(gps == levels(gps)[i]), which(gps == levels(gps)[j])]
      		S21 <- t(S12)
      		RV2[i,j] <- sum(diag(S12%*%S21))/sqrt(sum(diag(S11%*%S11))*sum(diag(S22%*%S22)))
      		diag(RV2) <- 0
    	}
 	 }
 	 RV2 <- sum(RV2)/(nmod/2*(nmod-1))
 	 
 	 # compute difference
 	 RV.d = RV1 - RV2
	
	# permute groups and compute RV1 and RV2
  	P.val <- 0
  	RV.val <- c()
  	for(ii in 1:iter) {
    	groups.r <- sample(groups)
    	if (length(dim(A)) == 3) {
			A1 <- A[,,which(groups.r == levels(groups.r)[1])]
			A1 <- gpagen(A1, ShowPlot=FALSE)[[1]]
			A2 <- A[,,which(groups.r == levels(groups.r)[2])]
			A2 <- gpagen(A2, ShowPlot=FALSE)[[1]]
			x1 <- two.d.array(A1)
			x2 <- two.d.array(A2)}
  		if (length(dim(A))==2){
  			x1 <- A[which(groups == levels(groups)[1]),]
        	x2 <- A[which(groups == levels(groups)[2]),]}
			S1 <- cov(x1)
  			RV1.r <- array(0, dim=c(nmod, nmod))
 			for (i in 1:(nmod-1)){
    			for (j in 2: nmod){
     				S11 <- S1[which(gps == levels(gps)[i]), which(gps == levels(gps)[i])]
      				S22 <- S1[which(gps == levels(gps)[j]), which(gps == levels(gps)[j])]
      				S12 <- S1[which(gps == levels(gps)[i]), which(gps == levels(gps)[j])]
      				S21 <- t(S12)
      				RV1.r[i,j] <- sum(diag(S12%*%S21))/sqrt(sum(diag(S11%*%S11))*sum(diag(S22%*%S22)))
      				diag(RV1.r) <- 0
    			}
 	 		}
 	 		RV1.r <- sum(RV1.r)/(nmod/2*(nmod-1))
			S2 <- cov(x2)
  			RV2.r <- array(0, dim=c(nmod, nmod))
 			for (i in 1:(nmod-1)){
    			for (j in 2:nmod){
     				S11 <- S2[which(gps == levels(gps)[i]), which(gps == levels(gps)[i])]
      				S22 <- S2[which(gps == levels(gps)[j]), which(gps == levels(gps)[j])]
      				S12 <- S2[which(gps == levels(gps)[i]), which(gps == levels(gps)[j])]
      				S21 <- t(S12)
      				RV2.r[i,j] <- sum(diag(S12%*%S21))/sqrt(sum(diag(S11%*%S11))*sum(diag(S22%*%S22)))
      				diag(RV2.r) <- 0
    			}
 	 		}
 	 		RV2.r <- sum(RV2.r)/(nmod/2*(nmod-1))
    		RV.d.r <- RV1.r - RV2.r
    		RV.val[ii] <- RV.d.r
    		if (RV.d.r < RV.d) {P.val <- P.val + 1}
    	}
    	
  RV.val[iter+1] <- RV.d
  P.val <- P.val/(iter+1)
  if (ShowPlot == TRUE) {
  	hist(RV.val,30,freq=TRUE,col="gray",xlab="RV1 - RV2", main=paste(levels(groups)[1], "vs.", levels(groups)[2]))
  	arrows(RV.d, 50, RV.d, 5, length=0.1, lwd=2)
  }
  
  return(data.frame(RV1, RV2, P.val))

}


