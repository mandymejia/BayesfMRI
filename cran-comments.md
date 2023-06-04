## New submission

This is the first submission of `BayesfMRI` to CRAN.

## Previous submission results

  Possibly misspelled words in DESCRIPTION:
    INLA (35:44, 36:75, 37:35, 38:5)
    dep (38:56)
    getOption (37:50)
    repos (37:42, 37:61)

We've removed the statement about INLA to avoid this warning.

  Suggests or Enhances not in mainstream repositories:
    INLA

Our package does require INLA for certain functionality. We've added the `Additional_repositories` statement for INLA to the DESCRIPTION.

  The Description field contains
    INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)`.
  Please enclose URLs in angle brackets (<...>).

We've removed the statement about INLA to avoid this warning.

## Test environments

* Windows x86_64-w64-mingw32/x64, R 4.2.2
* Mac aarch64-apple-darwin20 (64-bit), R 4.2.3

## R CMD check results

❯ checking compilation flags in Makevars ... WARNING
  Non-portable flags in variable 'PKG_CPPFLAGS':
    -O3
  Non-portable flags in variable 'PKG_CXXFLAGS':
    -O3

0 errors ✔ | 1 warning ✖ | 0 notes ✔

I'm sorry, I don't know what these mean. I can't find a "-03" flag in our package?

## Tests

Passes the roxygen `examples`, and the tests that require data (results not included in package) appear correct.
