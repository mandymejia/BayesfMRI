## Test environments

* Mac aarch64-apple-darwin20 (64-bit), R 4.4.0

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Tests

Passes the roxygen `examples`, and the tests that require data (results not included in package) appear correct.

## Previous submission

Suggests or Enhances not in mainstream repositories:
    INLA
  Availability using Additional_repositories specification:
    INLA   yes   https://inla.r-inla-download.org/R/testing

  Package has a VignetteBuilder field but no prebuilt vignette index.

* "VignetteBuilder: knitr" has been dropped from the DESCRIPTION for now, while there is not a vignette.

Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-x86_64
Check: S3 generic/method consistency, Result: NOTE
  Mismatches for apparent methods not registered:
  is.matrix:
    function(x)
  is.matrix.or.df:
    function(q)
  See section 'Registering S3 methods' in the 'Writing R Extensions'
  manual.

* This function has been renamed to `is_matrix_or_df`.

Flavor: r-devel-linux-x86_64-debian-gcc, r-devel-windows-x86_64
Check: Rd cross-references, Result: NOTE
  Found the following Rd file(s) with Rd \link{} targets missing package
  anchors:
    activations.posterior.Rd: excursions.inla
  Please provide package anchors for all Rd \link{} targets not in the
  package itself and the base packages.

* The link has been removed.