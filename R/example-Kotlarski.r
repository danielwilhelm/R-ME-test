set.seed(1090)
library(METests)

# generate true explanatory variable
Xstar <- runif(100)

# generate repeated measurements X and Z
X <- Xstar + rnorm(100, 0, 0.3)
Z <- Xstar + rnorm(100, 0, 0.3)

# generate outcome
Y <- Xstar^2+Xstar/2 + rnorm(100, 0, 0.5)	

# perform the test for measurement error using cross-validated bandwidth
DMTest(Y, X, Z)