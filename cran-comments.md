## Initial submission






## Test environments

* macOS 10.15.7 (local), R 4.0.3


* ubuntu 16.04 (on travis-ci), R-devel, R 4.0.0, R 3.6.3
* macOS 10.13.6 (on travis-ci), R 4.0.1, R 3.6.3
* Windows (win-builder), R-devel, R 4.0.0, R 3.6.3
* Rhub
    - Oracle Solaris 10, R-patched, 32-bit
    - Fedora Linux, R-devel, clang (with valgrind)



## R CMD check results


There were no ERRORs or WARNINGs.


There were 2 NOTEs:

```
checking installed package size ... NOTE
    installed size is  5.5Mb
    sub-directories of 1Mb or more:
      libs   5.0Mb
```

The package makes extensive use of compiled code for improved performance
of simulations and for the model code via `rstan`.


```
GNU make is a SystemRequirements.
```

GNU make syntax is required to for packages using `rstan`.




## Downstream dependencies

There are currently no downstream dependencies for this package
