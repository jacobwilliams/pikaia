                                                       20 January, 1996

UDPATED: 2025, the data files are now JSON.

This directory contains driver codes, fitness functions and synthetic
datasets for the examples discussed in Sect. 5 of the PIKAIA User's
guide.

Note that the drivers assume that the required files for synthetic
datasets reside in the test directory.
All files for the synthetic datasets were originally FORTRAN unformatted
IEEE single-precision (big_endian).
They were converted to JSON for better portability in 2025.
The drivers include an initialization
subroutine (finit) that reads in these data.

The following table lists the drivers and fitness functions, together
with the corresponding Section in the User's guide

----------------------------------------------------------------------
Section   Driver    Fitness function   Dataset
----------------------------------------------------------------------
5.2       xkp1a.f       fit1a.f       syndat1.json    [Figure 5.4]
5.3       xkp1b.f       fit1b.f       syndat1.json    [Figure 5.4]
5.4       xkp2.f        fit2.f        syndat2.json    [Figure 5.9]
5.5       xkp3.f        fit3.f        syndat3.json    [Figure 5.10]
----------------------------------------------------------------------

The driver and fitness function for the 2-D landscape problem
of Section 5.1 are include within the pikaia.f file in the
main pikaia directory.
