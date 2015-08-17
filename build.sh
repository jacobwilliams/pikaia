#!/bin/bash

#
#  NAME
#    build.sh
#
#  DESCRIPTION
#    Build the pikaia library and unit test.
#
#  USAGE
#    build.sh         : build using gfortran
#    build.sh -ifort  : build using ifort
#
#  REQUIRES
#    FoBiS.py : https://github.com/szaghi/FoBiS
#    FORD     : https://github.com/cmacmackin/ford
#
#  AUTHOR
#    Jacob Williams : 3/8/2015 (based on the one from json-fortran)
#

set -e

DOCDIR='./documentation/'       # build directory for documentation
SRCDIR='./src/'                 # library source directory
TESTDIR='./src/tests/'          # unit test source directory
BINDIR='./bin/'                 # build directory for unit tests
LIBDIR='./lib/'                 # build directory for library
MODCODE='pikaia_module.f90'     # module file name
LIBOUT='libpikaia.a'            # name of library

if [ "$1" == "-ifort" ]; then
    # Intel compiler

    FCOMPILER='Intel'
    # The following warning might be triggered by ifort unless explicitly silenced:
    # warning #7601: F2008 standard does not allow an internal procedure to be an actual argument procedure name. (R1214.4).
    # In the context of F2008 this is an erroneous warning.
    # See https://prd1idz.cps.intel.com/en-us/forums/topic/486629
    FCOMPILERFLAGS='-c -O2 -warn -stand f08 -diag-disable 7601 -traceback'
    #FCOMPILERFLAGS='-c -O2 -warn -traceback -stand f08 -assume protect_parens -assume buffered_io -check all'

else
    # GFortran

    FCOMPILER='gnu'
    FCOMPILERFLAGS='-c -O2 -fbacktrace -Wall -Wextra -Wno-maybe-uninitialized -pedantic -std=f2008'

fi

#build source using FoBiS:

if hash FoBiS.py 2>/dev/null; then

    #build the stand-alone library:
    echo ""
    echo "Building library..."
    FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${LIBDIR} -s ${SRCDIR} -dmod ./ -dobj ./ -t ${MODCODE} -o ${LIBOUT} -mklib static -colors

    #build the unit tests (uses the above library):
    if [[ $JF_SKIP_TESTS != [yY]* ]]; then
        echo ""
        echo "Building unit tests..."
        for TEST in "${TESTDIR%/}"/*.f90; do
        THIS_TEST=${TEST##*/}
        echo "Build ${THIS_TEST%.f90}"
        FoBiS.py build -compiler ${FCOMPILER} -cflags "${FCOMPILERFLAGS}" -dbld ${BINDIR} -s ${TESTDIR} -i ${LIBDIR} -libs ${LIBDIR}/${LIBOUT} -dmod ./ -dobj ./ -t ${THIS_TEST} -o ${THIS_TEST%.f90} -colors
        done
    else
        echo "Skip building the unit tests since \$JF_SKIP_TESTS has been set to 'true'."
    fi

else
    echo "FoBiS.py not found! Cannot build library. Install using: sudo pip install FoBiS.py"
fi

#build the documentation with FORD:

if hash ford 2>/dev/null; then
    echo "Building documentation..."
    ford ./pikaia.md
else
    echo "Ford not found! Cannot build documentation. Install using: sudo pip install ford"
fi
