#!/bin/sh

f90flags=
if [ $1 = "-t" ]; then
   # Note: gfortran-specific check flag.
   f90flags="$fflags -DTEST_MODE -fcheck=all"
fi

python -m numpy.f2py -c --f2cmap .f2py_f2cmap -m fort_zmtc --f90flags="$f90flags" shr_kind_mod.F90 zm_timeless_closure.F90
