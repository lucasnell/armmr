
## Test environments

* macOS 10.15.7 (local), R 4.0.3


* ubuntu 16.04 (travis-ci), R-devel, R 4.0.2, R 3.6.3
* macOS 10.13.6 (travis-ci), R 4.0.3, R 3.6.3

* Windows (win-builder), R-devel, R 4.0.0, R 3.6.3
* Oracle Solaris 10 (R-hub), R-patched, 32-bit
* Fedora Linux (R-hub), R-devel, clang (with valgrind)



## R CMD check results

There were no ERRORs or WARNINGs.

There were 2 NOTEs:

```
Maintainer: 'Lucas A. Nell <lucas@lucasnell.com>'

New submission
```
This is indeed a new submission for this package.


```
checking installed package size ... NOTE
    installed size is  5.5Mb
    sub-directories of 1Mb or more:
      libs   5.2Mb
```

The package makes extensive use of compiled code for improved performance
of simulations and for the model code via `rstan`.



## Downstream dependencies

There are currently no downstream dependencies for this package
