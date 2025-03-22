                                                       20 January, 1996

This directory contains driver codes, fitness functions and synthetic
datasets for the examples discussed in Sect. 5 of the PIKAIA User's
guide. To run any of the examples, one only needs to replace, in the
self-contained code pikaia.f you obtained from in the main pikaia
directory, the driver xpkaia and fitness function twod by the desired
driver and fitness function taken from this directory.

Note that the drivers assume that the required files for synthetic
datasets reside in the directory where the driver will be run.
All files for the synthetic datasets are FORTRAN unformatted
IEEE single-precision. The drivers include an initialization
subroutine (finit) that reads in these data.

The following table lists the drivers and fitness functions, together
with the corresponding Section in the User's guide

----------------------------------------------------------------------
Section   Driver    Fitness function   Dataset
----------------------------------------------------------------------
5.2       xkp1a.f       fit1a.f       syndat1.f    [Figure 5.4]
5.3       xkp1b.f       fit1b.f       syndat1.f    [Figure 5.4]
5.4       xkp2.f        fit2.f        syndat2.f    [Figure 5.9]
5.5       xkp3.f        fit3.f        syndat3.f    [Figure 5.10]
----------------------------------------------------------------------

The driver and fitness function for the 2-D landscape problem
of Section 5.1 are include within the pikaia.f file in the
main pikaia directory.
