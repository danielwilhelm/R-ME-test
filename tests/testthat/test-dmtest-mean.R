context("DM Test for conditional mean independence")

test_that("return value is of correct class and size", {

	set.seed(100)
	n <- 100
	X0 <- rnorm(n); Z0 <- rnorm(n); Y0 <- X0 + rnorm(n)

	for (stat in c("CvM", "KS")) {
		for (a in c(NA,0.5)) {
			res <- DMTest(Y0, X0, Z0, a=a, stat=stat)
			
			expect_type(res$teststat, "double")
			expect_type(res$cv, "double")
			expect_type(res$rej, "logical")
			expect_type(res$a, "double")
			expect_type(res$pval, "double")
			expect_lt(res$teststat, 1)
			expect_gte(res$teststat, 0)
			expect_false(res$rej)
			expect_gt(res$pval, 0.1)

		}
	}
})

