#!/bin/sh

f90flags="-lzmtc"
only_list="cape_consumption_rate_ffi"
if [ $1 = "-t" ]; then
    # Note: gfortran-specific check flag.
    f90flags="$fflags -DTEST_MODE -fcheck=all"
    only_list="weight_ffi weight_1d_ffi cape_consumption_ongoing_ffi cape_consumption_starting_ffi cape_consumption_ending_ffi end_time_frac_ffi $only_list"
fi

gfortran -shared -o libzmtc.so -fPIC $f90flags shr_kind_mod.F90 zm_timeless_closure.F90
export LDFLAGS="-Wl,-rpath=."
python -m numpy.f2py -c --f2cmap .f2py_f2cmap -m fort_zmtc --f90flags="$f90flags" -L. -lzmtc zm_timeless_closure_ffi.F90 only: $only_list
