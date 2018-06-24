# METests

The `R` package `METests` implements tests for various hypotheses about measurement error as proposed in Wilhelm (2017). The R package contains help files describing the various commands, their syntax and gives examples.

## Installation

1. Install the package `devtools` if it isn't already:

```R
install.packages("devtools")
```

2. Load the package `devtools`:

```R
library("devtools")
```

3. Install the package `METests`:

```R
install_github("danielwilhelm/R-ME-test")
```

## Example

Here is an example of how to test the null hypothesis that there is no measurement error in an explanatory variable X, using a second measurement Z and an outcome Y that depends on the true explanatory variable:

```R
rm(list = ls(all = TRUE))
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
```


# Reference
[Wilhelm, D. (2017), "Testing for the Presence of Measurement Error", Working Paper available soon](http://www.ucl.ac.uk/~uctpdwi/)