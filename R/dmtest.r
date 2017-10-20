#' compute Delgado and Manteiga (2001) test statistic
#'
#' compute Delgado and Manteiga (2001) test statistic for testing the hypothesis H0: E[Y|X,Z] = E[Y|X]
#'	
#' @param Y n-dim. vector containing the observations on the outcome
#' @param X matrix with n rows containing the observations on the variable X
#' @param Z matrix with n rows containing the observations on the variable Z
#' @param a vector of bandwidths, of the same dimension as there are columns in X, if unspecified, then the bandwidths are determined by cross-validation from nonparametric regression of Y on X
#' @param ckertype character string denoting the kernel function to be used, as in np package (default: "gaussian")
#' @param stat character string denoting the type of test statistic to be computed: Cramer-von-Mises ("CvM", default) or Kolmogorov-Smirnov ("KS")
#' @return (teststat, epsilonhat, Yhat) value of the test statistic, the residuals epsilonhat, and predicted values Yhat
#' @keywords significance test Delgado Manteiga measurement error
#' @export
#' @examples
#' computeTStat(Y, X, Z, a=NA, ckertype="gaussian", stat="CvM")
computeTStat <- function(Y, X, Z, a=NA, ckertype="gaussian", stat="CvM") {
	
	n <- length(Y)
	if (is.matrix(X) | is.data.frame(X)) dX <- ncol(X) else dX <- 1
	if (is.matrix(Z) | is.data.frame(Z)) dZ <- ncol(X) else dZ <- 1

	# nonparametric regression of Y on X
	if (is.na(a) | length(a)!=dX) { a <- np::npregbw(Y~X); print(a) }
	fhat <- c(np::npksum(txdat=X, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / n
	Yhat <- c(np::npksum(txdat=X, tydat=Y, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / fhat / n
	epsilonhat <- Y-Yhat

	# test statistic
	Tn <- function(x,z) mean(epsilonhat*fhat* (X<=x)*(Z<=z))
	Tn.vals <- c(apply(cbind(X,Z), 1, function(x) Tn(x[1:dX],x[(dX+1):(dX+dZ)])))
	teststat <- switch(stat,
		"CvM" = sum(Tn.vals^2),
		"KS"  = max(abs(sqrt(n)*Tn.vals)))

	return(list(teststat=teststat, epsilonhat=epsilonhat, Yhat=Yhat))
}


#' perform the Delgado and Manteiga (2001) test
#'
#' perform the Delgado and Manteiga (2001) test of the hypothesis H0: E[Y|X,Z] = E[Y|X]
#'	
#' @param Y n-dim. vector containing the observations on the outcome
#' @param X matrix with n rows containing the observations on the variable X
#' @param Z matrix with n rows containing the observations on the variable Z
#' @param size scalar between 0 and 1, denoting the nominal size of the test (default: 0.05)
#' @param B integer denoting the number of bootstrap samples to be used (default: 100)
#' @param a vector of bandwidths, of the same dimension as there are columns in X, if unspecified, then the bandwidths are determined by cross-validation from nonparametric regression of Y on X
#' @param ckertype character string denoting the kernel function to be used, as in np package (default: "gaussian")
#' @param stat character string denoting the type of test statistic to be computed: Cramer-von-Mises ("CvM", default) or Kolmogorov-Smirnov ("KS")
#' @return a list containing the following elements: 'teststat' value of the test statistic, 'cv' bootstrap critical value, 'rej' a 1-0 indicator for whether the test rejects or not, 'pval' p-value
#' @keywords significance test Delgado Manteiga measurement error
#' @export
#' @examples
#' DMTest(Y, X, Z, size=0.05, B=100, a=NA, ckertype="gaussian", stat="CvM")
DMTest <- function(Y, X, Z, size=0.05, B=100, a=NA, ckertype="gaussian", stat="CvM") {

	# compute test statistic
	res <- computeTStat(Y, X, Z, a, ckertype, stat)
	teststat <- res$teststat
	epsilonhat <- res$epsilonhat
	Yhat <- res$Yhat

	# bootstrap critical value
	teststatb <- rep(0,B)
	for (b in 1:B) {
		V <- rMammen(length(Y))
		Yhatstar <- Yhat+epsilonhat*V
		teststatb[b] <- computeTStat(Yhatstar, X, Z, a, ckertype, stat)$teststat
	}
	cv <- quantile(teststatb, 1-size)

	Fn <- ecdf(teststatb)
	pval <- 1-Fn(teststat)

	return(list(teststat=teststat, cv=cv, rej=teststat>cv, pval=pval))
}

#' draw a random sample from Mammen's two-point distribution
#'	
#' @param n sample size
#' @return V random sample from Mammen's two-point distribution
#' @keywords random sample generator Mammen
#' @export
#' @examples
#' rMammen(n)
rMammen <- function(n) {
	V <- rbinom(n, 1, prob=(sqrt(5)-1)/(2*sqrt(5)))
	V[V==0] <- (1-sqrt(5))/2
	V[V==1] <- (1+sqrt(5))/2
	return(V)
}

