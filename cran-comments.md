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
