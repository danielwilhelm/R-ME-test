#' compute Delgado and Manteiga (2001) test statistic
#'
#' compute Delgado and Manteiga (2001) test statistic for testing the hypothesis H0: E[Y|X,Z] = E[Y|X]
#'	
#' @param Y n-dim. vector containing the observations on the outcome
#' @param X matrix with n rows containing the observations on the scalar or vector X
#' @param Z matrix with n rows containing the observations on the scalar or vector Z
#' @param a vector of bandwidths, of the same dimension as there are columns in X, if unspecified, then the bandwidths are determined by cross-validation from nonparametric regression of Y on X
#' @param ckertype character string denoting the kernel function to be used, as in np package (default: "gaussian")
#' @param stat character string denoting the type of test statistic to be computed: Cramer-von-Mises ("CvM", default) or Kolmogorov-Smirnov ("KS")
#' @return a list containing the following elements: 'teststat' value of the test statistic, 'epsilonhat' the residuals from a nonparametric regression from Y on X, 'Yhat' the predicted values from a nonparametric regression of Y on X, 'a' the bandwidth(s)
#' @keywords significance test Delgado Manteiga measurement error
#' @export
#' @examples
#' Y <- rnorm(100)
#' X <- rnorm(100)
#' Z <- rnorm(100)
#' computeDMStat(Y, X, Z, a=NA, ckertype="gaussian", stat="CvM")
computeDMStat <- function(Y, X, Z, a=NA, ckertype="gaussian", stat="CvM") {

	HD <- dim(as.matrix(X))[2] > 1 | dim(as.matrix(Z))[2] > 1

	if (HD) {
		return(computeDMStat_HD(Y, X, Z, a=a, ckertype=ckertype, stat=stat))
	} else {
		K <- getKernel(ckertype)
		return(computeDMStat_1D(Y, c(X), c(Z), a=a, K=K, stat=stat))
	}
}


computeDMStat_HD <- function(Y, X, Z, a=NA, ckertype="gaussian", stat="CvM") {
	
	n <- length(Y)
	if (is.matrix(X) | is.data.frame(X)) dX <- ncol(X) else { stopifnot(is.numeric(X)); dX <- 1; }
	if (is.matrix(Z) | is.data.frame(Z)) dZ <- ncol(Z) else { stopifnot(is.numeric(Z)); dZ <- 1; }

	# nonparametric regression of Y on X
	if (any(is.na(a)) | length(a)!=dX) a <- (np::npregbw(xdat=X, ydat=Y))$bw
	fhat <- c(np::npksum(txdat=X, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / n
	Yhat <- c(np::npksum(txdat=X, tydat=Y, leave.one.out=FALSE, bandwidth.divide=TRUE, bws=a, ckertype=ckertype)$ksum) / fhat / n
	epsilonhat <- Y-Yhat

	# test statistic
	Tn <- function(x,z) {
		X <- as.matrix(X); Z <- as.matrix(Z)
		stopifnot(ncol(X)==length(x) & ncol(Z)==length(z))
		if (ncol(X)>1) x <- matrix(rep(x,n),nrow=n,byrow=TRUE)
		if (ncol(Z)>1) z <- matrix(rep(z,n),nrow=n,byrow=TRUE)
		ind <- apply((X<=x), 1, prod) * apply((Z<=z), 1, prod)
		return(mean(epsilonhat*fhat*ind))
	}
	Tn.vals <- c(apply(cbind(X,Z), 1, function(x) Tn(x[1:dX],x[(dX+1):(dX+dZ)])))
	teststat <- switch(stat,
		"CvM" = sum(Tn.vals^2),
		"KS"  = max(abs(sqrt(n)*Tn.vals)))

	return(list(teststat=teststat, epsilonhat=epsilonhat, Yhat=Yhat, a=a))
}

computeDMStat_1D <- function(Y, X, Z, a=NA, K=dnorm, stat="CvM") {
	n <- length(Y)

	# nonparametric regression of Y on X
	if (any(is.na(a)) | length(a)!=1) a <- (np::npregbw(xdat=X, ydat=Y))$bw
	Kij <- K(outer(X,X,"-")/a)/a
	fhat <- colMeans(Kij)
	Y.mat <- matrix(rep(Y,n), nrow=n, byrow=FALSE)
	Yhat <- colMeans(Y.mat*Kij) / fhat
	epsilonhat <- Y-Yhat

	# put into matrix form
	epsilonhat.mat <- matrix(rep(epsilonhat,n), nrow=n, byrow=FALSE)
	fhat.mat <- matrix(rep(fhat,n), nrow=n, byrow=FALSE)

	# compute test statistic	
	ind <- outer(X, X, "<=")*outer(Z, Z, "<=")
	Tn.vals <- colMeans(epsilonhat.mat*ind*fhat)

	teststat <- switch(stat,
		"CvM" = sum(Tn.vals^2),
		"KS"  = max(abs(sqrt(n)*Tn.vals)))

	return(list(teststat=teststat, epsilonhat=epsilonhat, Yhat=Yhat, a=a))
}

computeDMStat_1D_ind <- function(Y, X, Z, a=NA, K=dnorm, stat="CvM") {
	n <- length(Y)

	# nonparametric regression of Y on X
	if (any(is.na(a)) | length(a)!=1) a <- (np::npregbw(xdat=X, ydat=Y))$bw
	Kij <- K(outer(X,X,"-")/a)/a
	fhat <- colMeans(Kij)
	
	compTn <- function(y,x,z) {
		Y.mat <- matrix(rep(Y<=y,n), nrow=n, byrow=FALSE)	
		Yhat <- colMeans(Y.mat*Kij) / fhat
		epsilonhat <- Y-Yhat
		return(list(Tn=mean(epsilonhat*(X<=x)*(Z<=z)*fhat),
			Yhat=Yhat, epsilonhat=epsilonhat))
	}
	
	Tn.vals <- mapply(compTn, y=Y, x=X, z=Z)

	teststat <- switch(stat,
		"CvM" = sum(Tn.vals^2),
		"KS"  = max(abs(sqrt(n)*Tn.vals)))

	return(list(teststat=teststat, epsilonhat=epsilonhat, Yhat=Yhat, a=a))
}



#' perform the Delgado and Manteiga (2001) test
#'
#' perform the Delgado and Manteiga (2001) test of the hypothesis of conditional mean independence (E[Y|X,Z] = E[Y|X]) or conditional independence (Y is independent of Z given X)
#'	
#' @param Y n-dim. vector containing the observations on the outcome
#' @param X matrix with n rows containing the observations on the scalar or vector X
#' @param Z matrix with n rows containing the observations on the scalar or vector Z
#' @param size scalar between 0 and 1, denoting the nominal size of the test (default: 0.05)
#' @param B integer denoting the number of bootstrap samples to be used (default: 100)
#' @param a vector of bandwidths, of the same dimension as there are columns in X, if unspecified, then the bandwidths are determined by cross-validation from nonparametric regression of Y on X
#' @param ckertype character string denoting the kernel function to be used, as in np package (default: "gaussian")
#' @param stat character string denoting the type of test statistic to be computed: Cramer-von-Mises ("CvM", default) or Kolmogorov-Smirnov ("KS")
#' @param indep logical; FALSE means that conditional mean independence is tested, otherwise conditional independence
#' @return a list containing the following elements: 'teststat' value of the test statistic, 'cv' bootstrap critical value, 'rej' a 1-0 indicator for whether the test rejects or not, 'pval' p-value, 'a' the bandwidth(s)
#' @keywords significance test Delgado Manteiga measurement error
#' @export
#' @examples
#' Y <- rnorm(100)
#' X <- rnorm(100)
#' Z <- rnorm(100)
#' DMTest(Y, X, Z, size=0.05, B=100, a=NA, ckertype="gaussian", stat="CvM")
DMTest <- function(Y, X, Z, size=0.05, B=100, a=NA, ckertype="gaussian", stat="CvM", indep=FALSE) {

	stopifnot(size<1 & size>0)
	HD <- dim(as.matrix(X))[2] > 1 | dim(as.matrix(Z))[2] > 1
	if (indep & HD) stop("Test for conditional independence currently only implemented for scalar X and Z.")

	if (HD) {
		compDMStat <- function(ybar, abar) computeDMStat_HD(ybar, X, Z, a=abar, ckertype=ckertype, stat=stat)
	} else {
		K <- getKernel(ckertype)
		if (indep) {
			compDMStat <- function(ybar, abar) computeDMStat_1D_ind(ybar, c(X), c(Z), a=abar, K=K, stat=stat)
		}
		compDMStat <- function(ybar, abar) computeDMStat_1D(ybar, c(X), c(Z), a=abar, K=K, stat=stat)
	}

	# compute test statistic
	res <- compDMStat(Y, a)
	teststat <- res$teststat
	epsilonhat <- res$epsilonhat
	Yhat <- res$Yhat
	ahat <- res$a

	# bootstrap critical value
	teststatb <- replicate(B, compDMStat(epsilonhat*rMammen(length(Y)), ahat)$teststat)
	cv <- quantile(teststatb, 1-size, na.rm=TRUE)

	Fn <- ecdf(teststatb)
	pval <- 1-Fn(teststat)

	return(list(teststat=teststat, cv=cv, rej=teststat>cv, pval=pval, a=ahat))
}

#' draw a random sample from Mammen's two-point distribution
#'	
#' @param n sample size
#' @return V random sample from Mammen's two-point distribution
#' @keywords random sample generator Mammen
#' @export
#' @examples
#' rMammen(100)
rMammen <- function(n) {
	V <- rbinom(n, 1, prob=(sqrt(5)-1)/(2*sqrt(5)))
	V[V==0] <- (1-sqrt(5))/2
	V[V==1] <- (1+sqrt(5))/2
	return(V)
}


getKernel <- function(name="gaussian") {
	return(switch(name,
			"gaussian" = dnorm,
			"BT" = k_BT,
			"TR" = k_TR,
			"PR" = k_PR,
			"TH" = k_TH,
			"QS" = k_QS))
}

# Bartlett Kernel
k_BT <- function(x, renorm=FALSE) {
	if (renorm) x <- 2/3*x
	return(ifelse(abs(x)<=1,1-abs(x),0))
}

# Truncated Kernel
k_TR <- function(x, renorm=FALSE) {
	if (renorm) x <- 2*x
	return(ifelse(abs(x)<=1,1,0))
}

# Parzen Kernel
k_PR <- function(x, renorm=FALSE) {
	if (renorm) x <- 0.539285*x
	return(ifelse(abs(x) >= 0 & abs(x)<=0.5, 1-6*x^2+6*abs(x)^3, ifelse(abs(x) <= 1 & abs(x)>=0.5, 2*(1-abs(x))^3 ,0)))
}

# Tukey-Hanning Kernel
k_TH <- function(x, renorm=FALSE) {
	if (renorm) x <- 3/4*x
	return(ifelse(abs(x)<=1,(1+cos(pi*x))/2,0))
}

# Quadratic Spectral Kernel
k_QS <- function(x) return(25/(12*pi^2*x^2)*(sin(6*pi*x/5)/(6*pi*x/5)-cos(6*pi*x/5)))
