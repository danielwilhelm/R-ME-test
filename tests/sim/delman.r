library(np)

delman <- function(Y, X, Z, a=NA, ckertype="gaussian", stat="CvM") {
	
	n <- length(Y)
	if (is.matrix(X) | is.data.frame(X)) dX <- ncol(X) else dX <- 1
	if (is.matrix(Z) | is.data.frame(Z)) dZ <- ncol(X) else dZ <- 1

	# nonparametric regression of Y on X
	if (is.na(a) | length(a)!=dX) { a <- npregbw(Y~X); print(a) }
	fhat <- c(npksum(txdat=X, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / n
	Yhat <- c(npksum(txdat=X, tydat=Y, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / fhat / n
	epsilonhat <- Y-Yhat

	# test statistic
	Tn <- function(x,z) mean(epsilonhat*fhat* (X<=x)*(Z<=z))
	Tn.vals <- c(apply(cbind(X,Z), 1, function(x) Tn(x[1:dX],x[(dX+1):(dX+dZ)])))
	teststat <- switch(stat,
		"CvM" = sum(Tn.vals^2),
		"KS"  = max(abs(sqrt(n)*Tn.vals)))

	return(list(teststat=teststat, epsilonhat=epsilonhat, Yhat=Yhat))
}

# delman.cdist <- function(Y, X, Z, a=NA, ckertype="gaussian", stat="KS") {
	
# 	n <- length(Y)
# 	if (is.matrix(X) | is.data.frame(X)) dX <- ncol(X) else dX <- 1
# 	if (is.matrix(Z) | is.data.frame(Z)) dZ <- ncol(X) else dZ <- 1

# 	ygrid <- seq(quantile(Y,prob=0.1), quantile(Y,prob=0.9), length=100)

# 	# nonparametric regression of 1{Y<y} on X
# 	indY <- (Y <=)
# 	if (is.na(a) | length(a)!=dX) { a <- npregbw(Y~X); print(a) }
# 	fhat <- c(npksum(txdat=X, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / n
# 	Yhat <- c(npksum(txdat=X, tydat=Y, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / fhat / n
# 	epsilonhat <- Y-Yhat

# 	# test statistic
# 	Tn <- function(x,z) mean(epsilonhat*fhat* (X<=x)*(Z<=z))
# 	Tn.vals <- c(apply(cbind(X,Z), 1, function(x) Tn(x[1:dX],x[(dX+1):(dX+dZ)])))
# 	teststat <- max(abs(sqrt(n)*Tn.vals))

# 	return(list(teststat=teststat, epsilonhat=epsilonhat))
# }


delman.test <- function(Y, X, Z, size=0.05, B=100, a=NA, ckertype="gaussian", stat="CvM") {

	# compute test statistic
	res <- delman(Y, X, Z, a, ckertype, stat)
	teststat <- res$teststat
	epsilonhat <- res$epsilonhat
	Yhat <- res$Yhat

	# bootstrap critical value
	teststatb <- rep(0,B)
	for (b in 1:B) {
		V <- rMammen(length(Y))
		Yhatstar <- Yhat+epsilonhat*V
		teststatb[b] <- delman(Yhatstar, X, Z, a, ckertype, stat)$teststat
	}
	cv <- quantile(teststatb, 1-size)

	return(list(teststat=teststat, cv=cv, rej=teststat>cv))
}

rMammen <- function(n) {
	V <- rbinom(n, 1, prob=(sqrt(5)-1)/(2*sqrt(5)))
	V[V==0] <- (1-sqrt(5))/2
	V[V==1] <- (1+sqrt(5))/2
	return(V)
}

ttest <- function(Y, X, Z, size=0.05) {
	res <- lm(Y~X+Z)
	tstat <- coef(summary(res))["Z", "t value"]
	return(list(teststat=abs(tstat), cv=qnorm(1-size/2), rej=abs(tstat)>qnorm(1-size/2)))
}