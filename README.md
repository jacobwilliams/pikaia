![pikaia](media/logo.png)
============

# pikaia [![GitHub release](https://img.shields.io/github/release/jacobwilliams/pikaia.svg?style=plastic)](https://github.com/jacobwilliams/pikaia/releases/latest)
Modern Fortran Edition of the Pikaia Genetic Algorithm
### Status

![Build Status](https://github.com/jacobwilliams/pikaia/actions/workflows/CI.yml/badge.svg)

### Overview

This is a refactoring of the PIKAIA unconstrained optimization code from the [High Altitude Observatory](http://www.hao.ucar.edu/modeling/pikaia/pikaia.php).  The original code is public domain and was written by Paul Charbonneau & Barry Knapp.  The new code differs from the old code in the following respects:
 * The original fixed-form source (FORTRAN 77) was converted to free-form source.
 * The code is now object-oriented Fortran 2003/2008.  All user interaction is now through the ```pikaia_class```.
 * All real variables are now double precision.
 * The original random number generator was replaced with MT19937-64 (64-bit Mersenne Twister).
 * There are various new options (e.g., a convergence window with a tolerance can be specified as a stopping condition, and the user can specify a subroutine for reporting iterations).
 * Mapping the variables to be between 0 and 1 now occurs internally, rather than requiring the user to do it.
 * Can now include an initial guess in the initial population.
 * Some OpenMP support has been added.

### Compiling

A [FoBiS](https://github.com/szaghi/FoBiS) configuration file (`pikaia.fobis`) is also provided that can also build the library and examples. Use the `mode` flag to indicate what to build. For example:

  * To build all the examples using gfortran: `FoBiS.py build -f pikaia.fobis -mode tests-gnu`
  * To build all the examples using ifort: `FoBiS.py build -f pikaia.fobis -mode tests-intel`
  * To build a static library using gfortran: `FoBiS.py build -f pikaia.fobis -mode static-gnu`
  * To build a static library using ifort: `FoBiS.py build -f pikaia.fobis -mode static-intel`

  The full set of modes are: `static-gnu`, `static-gnu-debug`, `static-intel`, `static-intel-debug`, `shared-gnu`, `shared-gnu-debug`, `shared-intel`, `shared-intel-debug`, `tests-gnu`, `tests-gnu-debug`, `tests-intel`, `tests-intel-debug`,
  `static-gnu-openmp`, `static-gnu-debug-openmp`, `static-intel-openmp`, `static-intel-debug-openmp`, `shared-gnu-openmp`, `shared-gnu-debug-openmp`, `shared-intel-openmp`, `shared-intel-debug-openmp`, `tests-gnu-openmp`, `tests-gnu-debug-openmp`, `tests-intel-openmp`, `tests-intel-debug-openmp`

  To generate the documentation using [ford](https://github.com/Fortran-FOSS-Programmers/ford), run: ```FoBis.py rule --execute makedoc -f pikaia.fobis```

### Examples

 * An example use of Pikaia can be found [here](http://degenerateconic.com/earth-mars-free-return/).

### Documentation

 * The API documentation for the current ```master``` branch can be found [here](https://jacobwilliams.github.io/pikaia/).  This is generated by processing the source files with [FORD](https://github.com/Fortran-FOSS-Programmers/ford).  Note that the shell script will also generate these files automatically in the ```doc``` folder, assuming you have FORD installed.
 * The original Pikaia documentation (for v1.2) can be found [here](http://www.hao.ucar.edu/modeling/pikaia/relnotes.ps).

### See also

 * PIKAIA description page: http://www.hao.ucar.edu/modeling/pikaia/pikaia.php
 * Original source code: http://download.hao.ucar.edu/archive/pikaia/

