import os

from fort_zmtc import zm_timeless_closure as zmtc

if os.environ["TEST_MODE"] == "TRUE":
    weight = zmtc.weight
    weight_1d = zmtc.weight_1d
    cape_consumption_ongoing = zmtc.cape_consumption_ongoing
