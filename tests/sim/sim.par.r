###########
## SETUP ##
###########

rm(list = ls())
SAVE <- TRUE
model <- "Kotlarski"
mfun <- "quad"
# setwd("~/work/Documents/projects/testing_ME/sim/trunk")
# setwd("/Users/uctpdwi/Dropbox/work/sync/Documents/projects/testing_ME/sim/trunk")
# setwd(paste("/home/uctpdwi/testing_ME/",model,sep=""))
# source("delman.r")
library(doMC)
library(foreach)

# data-generating process

	dgp <- function(nobs=1000, beta=1, gammaX=1, gammaZ=1, sd.eta=0.5, model="Kotlarski") {
		
		if (mfun=="sine") m <- function(x) return(sin(pi*x/2))
		if (mfun=="hump") {
			m <- function(x) {
				m1 <- sin(pi*x)
				m2 <- -sin(2*pi*x)/(2^2)
				m3 <- sin(3*pi*x)/(3^2)
				m4 <- -sin(4*pi*x)/(4^2)
				return(m1+m2+m3+m4)	
			}
		}
		if (mfun=="quad") m <- function(x) return(x^2+x/2)


		if (model=="DM") {
			X <- runif(nobs)
			Z <- runif(nobs)
			U <- rnorm(nobs)
			m <- function(x) return(1+sin(10*x))
			Y <- m(X) + beta*sin(gammaZ*Z) + U
		}

		if (model=="quad") {
			Xstar <- runif(nobs)
			# X <- rnorm(nobs)
			# Z <- rnorm(nobs)
			# Y <- sd.eta*Z + X + rnorm(nobs, 0, sd.err)
			D <- rbinom(nobs, 1, beta)
			X <- Xstar + D*rnorm(nobs, 0, sd.eta)
			Z <- -(Xstar-1)^2 + rnorm(nobs, 0, 0.2)
			Y <- m(Xstar) + rnorm(nobs, 0, 0.2)	
		}

		if (model=="Kotlarski") {
			Xstar <- runif(nobs)
			# X <- rnorm(nobs)
			# Z <- rnorm(nobs)
			# Y <- sd.eta*Z + X + rnorm(nobs, 0, sd.err)
			D <- rbinom(nobs, 1, beta)
			X <- Xstar + D*rnorm(nobs, 0, sd.eta)
			Z <- Xstar + rnorm(nobs, 0, 0.3)
			Y <- m(Xstar) + rnorm(nobs, 0, 0.5)	
		}

		if (model=="nonclass") {
			Xstar <- runif(nobs)
			# X <- rnorm(nobs)
			# Z <- rnorm(nobs)
			# Y <- sd.eta*Z + X + rnorm(nobs, 0, sd.err)
			D <- rbinom(nobs, 1, beta)
			X <- Xstar + D*rnorm(nobs, 0, sd.eta)*exp(-abs(Xstar-1/2))
			Z <- Xstar + rnorm(nobs, 0, 0.3)*exp(-abs(Xstar-1/2))
			Y <- m(Xstar) + rnorm(nobs, 0, 0.5)	
		}

		if (model=="nonclass2") {
			Xstar <- runif(nobs)
			# X <- rnorm(nobs)
			# Z <- rnorm(nobs)
			# Y <- sd.eta*Z + X + rnorm(nobs, 0, sd.err)
			D <- rbinom(nobs, 1, beta)
			X <- Xstar + D*rnorm(nobs, 0, sd.eta)*exp(-abs(Xstar-1/2))
			Z <- Xstar + rnorm(nobs, 0, 0.3)
			Y <- m(Xstar) + rnorm(nobs, 0, 0.5)	
		}

		if (model=="panel") {
			X1star <- runif(nobs)
			X2star <- 1.5*X1star + rnorm(nobs, 0, 0.5)
			# X <- rnorm(nobs)
			# Z <- rnorm(nobs)
			# Y <- sd.eta*Z + X + rnorm(nobs, 0, sd.err)
			D <- rbinom(nobs, 1, beta)
			X <- X2star + D*rnorm(nobs, 0, sd.eta)
			Z <- X1star + rnorm(nobs, 0, 0.5)
			Y <- m(X2star) + rnorm(nobs, 0, 0.5)	
		}
		
		return(data.frame(Y=Y,X=X,Z=Z))
	}


# simulate power
	
	testpower <- function(nrep=100, size=0.05, test=c("DM opt", "DM-", "DM"), a=NA, ...) {
		rej <- matrix(rep(NA, length(test) * nrep), ncol = length(test))
		colnames(rej) <- test

		# find optimal bandwidth
		if (is.na(a) & ("DM opt" %in% test)) {
			a <- rep(0,5)
			for (cv in 1:5) {
				dat <- dgp(...)
				Y <- dat$Y^(1+(model=="hump"))
				X <- dat$X
				a[cv] <- npregbw(Y~X)$bw
			}
			aopt <- median(a)
		}
		

		# MC loop
		for (i in 1:nrep) {
			dat <- dgp(...)
			n <- length(dat$Y)
			compute_rej <- function(test) {
				test <- match.arg(test, c("DM opt", "DM-", "DM", "DM+", "t", "DMind"))
				res <- switch(test,
					"DM opt" 		= DMTest(Y=dat$Y^(1+(model=="hump")), X=dat$X, Z=dat$Z, size=size, a=aopt, ckertype="gaussian"),
					"DM-" 		= DMTest(Y=dat$Y^(1+(model=="hump")), X=dat$X, Z=dat$Z, size=size, a=0.1*n^(-1/3), ckertype="gaussian"),
					"DM" 	= DMTest(Y=dat$Y^(1+(model=="hump")), X=dat$X, Z=dat$Z, size=size, a=0.2*n^(-1/3), ckertype="gaussian"),
					"DM+" 	= DMTest(Y=dat$Y^(1+(model=="hump")), X=dat$X, Z=dat$Z, size=size, a=0.5*n^(-1/3), ckertype="gaussian"),
					"DMind"	= DMTest(Y=dat$Y^(1+(model=="hump")), X=dat$X, Z=dat$Z, size=size, a=0.2*n^(-1/3), ckertype="gaussian", indep=TRUE),
					"t" = ttest(Y=dat$Y, X=dat$X, Z=dat$Z, size=size))
				return(res$rej)
			}
			rej[i,] <- sapply(test, compute_rej)
		}
		return(colMeans(rej))
	}


# simulation loop

	simulation <- function(nobs=100, sd.eta=seq(0,2,length=10), gammaZ=c(0.1,1), beta=c(0.1,1), test=c("DM opt", "DM-", "DM"), a=NA, ...) {
		
		# simulate power of the tests
		print("simulate the power of the tests ...")
		prs <- expand.grid(nobs=nobs,sd.eta=sd.eta, gammaZ=gammaZ, beta=beta)
		nprs <- nrow(prs)
		ntest <- length(test)


		registerDoMC(detectCores()-(1-SAVE))
		pow <- matrix(rep(NA, ntest * nprs), ncol=ntest)
		# pb <- txtProgressBar(min = 1, max = nprs, style = 3)
		pow <- foreach (i=1:nprs, .combine=rbind) %dopar% {
			# print(setTxtProgressBar(pb, i))
			sink("log.txt", append=TRUE)
			print(paste(i, "/", nprs))
			testpower(test=test, nobs=prs$nobs[i], beta=prs$beta[i], gammaZ=prs$gammaZ[i], sd.eta=prs$sd.eta[i], a=a, ...)
		}
		# close(pb)
		# dev.off()


		# construct return value for power simulations
		rvalpow <- data.frame()
		for(i in 1:ntest) rvalpow <- rbind(rvalpow, prs)
		rvalpow$test <- gl(ntest, nprs, labels = test)
		rvalpow$power <- as.vector(pow)
		rvalpow$sd.eta <- rep(factor(prs$sd.eta), ntest)
		rvalpow$gammaZ <- rep(factor(prs$gammaZ), ntest)
		rvalpow$beta <- rep(factor(prs$beta), ntest)
		rvalpow$nobs <- rep(factor(prs$nobs), ntest)
		
		return(rvalpow)
	}



####################
## RUN & EVALUATE ##
####################

## run simulation
set.seed(1090)
starting_time <- proc.time()
sc_sim <- simulation(nrep=100, nobs=500, beta=seq(0,1,length=5), sd.eta=c(0.2,0.5,1), gammaZ=1, test=c("DMind"), a=NA, model=model)
# sc_sim <- simulation(nrep=1000, nobs=c(200,500), beta=seq(0,1,length=5), sd.eta=c(0.2,0.5,1), gammaZ=1, test=c("t", "DM-", "DM", "DM+"), a=NA, model=model)
print((proc.time()-starting_time)[1])

## save work space
if (SAVE) save.image(paste("workspace_",model,".Rdata",sep=""))

## numerical summary
tab <- xtabs(power ~ sd.eta + beta + nobs + test, data = sc_sim)
restable <- ftable(tab, row.vars = c("nobs", "sd.eta", "test"), col.vars = "beta")
print(restable)
if (SAVE) { library(xtable); source("ftable2data.frame.r"); print.xtable(ftable2data.frame(restable), type="latex", file=paste("table_",model,".tex",sep="")); }

## graphical summary
library("lattice")
bwtheme <- standard.theme("pdf", color=FALSE)
if (SAVE) pdf(file=paste("power_",model,".pdf", sep=""),width=10,height=8) else quartz(width=14,height=12)
xyplot(power ~ beta |  nobs + sd.eta, groups = ~ test, data = sc_sim, type = "b", ylim=c(-0.1,1.1), ylab="", main="rejection probability", auto.key=list(space="right", cex.title=1, lines=TRUE, points=TRUE), par.settings=bwtheme)
if (SAVE) dev.off()
 



