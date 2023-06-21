## New submission

This is the first submission of `BayesfMRI` to CRAN.

## Previous submission results (0.3.1)

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

## Previous submission results (0.3.2)

  Non-portable flags in variable 'PKG_CPPFLAGS':
    -O3
  Non-portable flags in variable 'PKG_CXXFLAGS':
    -O3

We've removed `PKG_CPPFLAGS` and `PKG_CXXFLAGS` from `src/Makevars`.

## Previous submission results (0.3.3)

  Please always explain all acronyms in the description text. -> fMRI

Done!

  If there are references describing the methods in your package, please
  add these in the description field of your DESCRIPTION file in the form
  authors (year) <doi:...>
  authors (year) <arXiv:...>
  authors (year, ISBN:...)
  or if those are not available: <https:...>
  with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for
  auto-linking. (If you want to add a title as well please put it in
  quotes: "Title")

Done!

  Please add \value to .Rd files regarding exported methods and explain
  the functions results in the documentation. Please write about the
  structure of the output (class) and also what the output means. (If a
  function does not return a value, please document that too, e.g.
  \value{No return value, called for side effects} or similar)
  Missing Rd-tags:
        id_activations.Rd: \value
        plot.act_BayesGLM_cifti.Rd: \value
        plot.BayesGLM_cifti.Rd: \value
        plot.BayesGLM2_cifti.Rd: \value
        summary.act_BayesGLM.Rd: \value
        summary.act_BayesGLM_cifti.Rd: \value
        summary.BayesGLM.Rd: \value
        summary.BayesGLM_cifti.Rd: \value
        summary.BayesGLM2.Rd: \value
        summary.BayesGLM2_cifti.Rd: \value

Done!

  Please add small executable examples in your Rd-files to illustrate the
  use of the exported function but also enable automatic testing.

Added one more, but the other functions cannot be run as an example without
  providing several large data objects. 

  You write information messages to the console that cannot be easily
  suppressed.
  It is more R like to generate objects that can be used to extract the
  information a user is interested in, and then print() that object.
  Instead of print()/cat() rather use message()/warning() or
  if(verbose)cat(..) (or maybe stop()) if you really have to write text to
  the console. (except for print, summary, interactive functions)

Wrapped the calls to `cat()` with `if (verbose)`.

  Please ensure that you do not use more than 2 cores in your examples,
  vignettes, etc.

No example uses more than two cores, and we don't have vignettes at this point.

## Previous submission results (0.3.4)

  We still see no \value in:
  Missing Rd-tags:
        plot.act_BayesGLM_cifti.Rd: \value
        plot.BayesGLM_cifti.Rd: \value
        plot.BayesGLM2_cifti.Rd: \value

Fixed!

## Test environments

* Windows x86_64-w64-mingw32/x64, R 4.2.2
* Mac aarch64-apple-darwin20 (64-bit), R 4.2.3

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Tests

Passes the roxygen `examples`, and the tests that require data (results not included in package) appear correct.
