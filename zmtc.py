import os

from fort_zmtc import zm_timeless_closure_ffi as zmtc

if "TEST_MODE" in os.environ and os.environ["TEST_MODE"] == "TRUE":
    weight = zmtc.weight_ffi
    weight_1d = zmtc.weight_1d_ffi
    cape_consumption_ongoing = zmtc.cape_consumption_ongoing_ffi
    cape_consumption_starting = zmtc.cape_consumption_starting_ffi
    cape_consumption_ending = zmtc.cape_consumption_ending_ffi
    end_time_frac = zmtc.end_time_frac_ffi

cape_consumption_rate = zmtc.cape_consumption_rate_ffi
