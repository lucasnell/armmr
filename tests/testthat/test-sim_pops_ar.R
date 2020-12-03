context("test-sim_pops_ar.R")

test_that("simulating communities", {
    n_spp <- 3
    max_t <- 10
    n_locs <- 2
    X <- matrix(rlnorm(max_t * n_locs), max_t, n_locs)
    N0 <- matrix(log(10), n_spp, n_locs)
    b0 <- matrix(log(100), n_spp, n_locs)
    b1 <- matrix(0.1, n_spp, n_locs)
    rho <- matrix(0.2, n_spp, n_locs)
    vcv <- diag(n_spp)
    vcv[lower.tri(vcv)] <- vcv[upper.tri(vcv)] <- 0.1
    vcv <- replicate(n_locs, vcv, simplify = FALSE)
    obs <- rep(0.1, n_spp)
    N <- sim_pops_ar(X, N0, b0, b1, rho, vcv, obs)
    class(N)
    expect_is(N, "list")
    expect_true(all(sapply(N, inherits, what = "matrix")))
})
