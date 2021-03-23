test_that("Miscellaneous functions are working", {
  if (requireNamespace("neuRosim", quietly = TRUE)) {
    q <- seq(0, 30, by=1/100)
    testthat::expect_equal(BayesfMRI:::canonical_HRF(q), neuRosim::canonicalHRF(q))
  }
})
