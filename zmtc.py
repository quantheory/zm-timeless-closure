import os

from fort_zmtc import zm_timeless_closure as zmtc

if os.environ["TEST_MODE"] == "TRUE":
    weight = zmtc.weight
    weight_1d = zmtc.weight_1d
    cape_consumption_ongoing = zmtc.cape_consumption_ongoing
    cape_consumption_starting = zmtc.cape_consumption_starting
    cape_consumption_ending = zmtc.cape_consumption_ending
    end_time_frac = zmtc.end_time_frac

cape_consumption_rate = zmtc.cape_consumption_rate
