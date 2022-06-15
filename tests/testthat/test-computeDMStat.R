context("computing DM test statistic")

test_that("return value is of correct class and size", {

	set.seed(100)
	n <- 100
	X0 <- rnorm(n); Z0 <- rnorm(n); Y0 <- X0 + rnorm(n)

	for (stat in c("CvM", "KS")) {
		for (a in c(NA,0.5)) {
			res <- computeDMStat(Y0, X0, Z0, a=a, stat=stat)
			
			expect_type(res$teststat, "double")
			expect_type(res$Yhat, "double")
			expect_type(res$epsilonhat, "double")
			expect_type(res$a, "double")
			expect_true(res$teststat<=1 & res$teststat >= 0)
			expect_equal(length(res$epsilonhat), n)
			expect_equal(length(res$Yhat), n)
		}
	}
})

