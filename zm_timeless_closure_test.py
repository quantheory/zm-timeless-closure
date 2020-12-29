#!/usr/bin/env python3

import os
os.environ["TEST_MODE"] = "TRUE"

import unittest

import numpy as np

import zmtc

class VectorCase(unittest.TestCase):
    """
    A TestCase with added methods for testing vectors of floats.

    Methods
    -------
    assertAlmostEqualVec(first, second, msg=None, **kwargs)
    """

    def assertAlmostEqualVec(self, first, second, msg=None, **kwargs):
        """A vector version of assertAlmostEqual.

        The lengths of `first` and `second` are asserted to be equal, then each
        element of `first` is asserted to be almost equal to the corresponding
        element of `second`.

        Any keyword arguments other than `msg` are passed on to the underlying
        assertAlmostEqual call.

        Parameters
        ----------
        first : iterable object
            First vector to compare.
        second : iterable object
            Second vector to compare.
        msg : str, optional
            Message passed down to assert method calls.
        """
        if msg is None:
            msg = ""
        else:
            msg = "; {}".format(msg)
        self.assertEqual(len(first), len(second),
                         msg="arrays of unequal size{}".format(msg))
        for i in range(len(first)):
            self.assertAlmostEqual(first[i], second[i],
                                   msg="unequal at index {}{}".format(i, msg),
                                   **kwargs)


class TestWeight(VectorCase):

    def test_weight(self):
        """Scalar weight function uses the correct formula."""
        timestep = 0.1
        weight = (1. - np.exp(-timestep)) / timestep
        zmtc_weight = zmtc.weight(timestep)
        self.assertAlmostEqual(weight, zmtc_weight)

    def test_weight_1d(self):
        """Vector weight function uses the correct formula."""
        timesteps = np.array([0.1, 0.5, 1.0])
        weights = (1. - np.exp(-timesteps)) / timesteps
        zmtc_weights = zmtc.weight_1d(timesteps)
        self.assertAlmostEqualVec(weights, zmtc_weights)

class TestConvectionEndingFraction(VectorCase):

    def test_convection_ending_fraction(self):
        """Time to end of convection uses the correct formula."""
        tau = 2700.
        deltat = 1800.
        a_crit = 1.
        a_prev = np.array([1.3, 1.5, 2., 1.5])
        a_star = np.array([0.25, 1., 1.01, 1.4])
        f_fac = (a_prev - a_crit) / (a_prev - a_star)
        f = tau * np.log(1 + f_fac*deltat / tau) / deltat
        zmtc_f = zmtc.end_time_frac(tau, a_crit, deltat, a_prev, a_star)
        self.assertAlmostEqualVec(f, zmtc_f)


# Have less tolerance for error for d(CAPE)/dt, since with the units we use here
# the values will be on the order of 10^-3 rather than 10^0.
DCAPE_PLACES = 10

class TestCapeConsumptions(VectorCase):

    def test_ongoing_cape_consumption(self):
        """Ongoing CAPE consumption function uses the correct formula."""
        tau = 2700.
        deltat = 1800.
        weight = zmtc.weight(deltat/tau)
        a_prev = np.array([1.25, 2., 1.5])
        a_star = np.array([1.5, 2., 1.4])
        a_crit = 1.
        cape_consumed = (a_prev - a_crit) * weight / tau + \
                        (a_star - a_prev) * (1 - weight) / deltat
        zmtc_cape_consumed = zmtc.cape_consumption_ongoing(tau, a_crit, deltat,
                                                           weight, a_prev,
                                                           a_star)
        self.assertAlmostEqualVec(cape_consumed, zmtc_cape_consumed,
                                  places=DCAPE_PLACES)

    def test_starting_cape_consumption(self):
        """Starting CAPE consumption function uses the correct formula."""
        tau = 2700.
        deltat = 1800.
        a_crit = 1.
        a_prev = np.array([0.25, 1., 0.5])
        a_star = np.array([1.3, 1.5, 2.])
        f = (a_crit - a_prev) / (a_star - a_prev)
        weight = zmtc.weight_1d((1-f)*deltat/tau)
        cape_consumed = (a_star - a_crit) * (1. - weight) / deltat
        zmtc_cape_consumed = zmtc.cape_consumption_starting(tau, a_crit, deltat,
                                                            a_prev, a_star)
        self.assertAlmostEqualVec(cape_consumed, zmtc_cape_consumed,
                                  places=DCAPE_PLACES)

    def test_ending_cape_consumption(self):
        """Ending CAPE consumption function uses the correct formula."""
        tau = 2700.
        deltat = 1800.
        a_crit = 1.
        a_prev = np.array([1.3, 1.5, 2.])
        a_star = np.array([0.25, 1., 1.01])
        f_fac = (a_prev - a_crit) / (a_prev - a_star)
        f = tau * np.log(1 + f_fac*deltat / tau) / deltat
        assert all(f < 1.), "incorrect input for ending CAPE consumption"
        cape_consumed = (f * a_star + (1. - f) * a_prev - a_crit) / deltat
        zmtc_cape_consumed = zmtc.cape_consumption_ending(a_crit, deltat, f,
                                                          a_prev, a_star)
        self.assertAlmostEqualVec(cape_consumed, zmtc_cape_consumed,
                                  places=DCAPE_PLACES)


if __name__ == "__main__":
    unittest.main()
