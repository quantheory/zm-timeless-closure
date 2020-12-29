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

class TestCapeConsumptions(VectorCase):

    def test_ongoing_cape_consumption(self):
        """Ongoing CAPE consumption function uses the correct formula."""
        tau = 1800
        deltat = 1800
        weight = zmtc.weight(deltat/tau)
        a_prev = np.array([1.25, 2., 1.5])
        a_star = np.array([1.5, 3., 2.])
        a_crit = 1.
        cape_consumed = (a_prev - a_crit) * weight / tau + \
                        (a_star - a_prev) * (1 - weight) / deltat
        zmtc_cape_consumed = zmtc.cape_consumption_ongoing(tau, a_crit, deltat,
                                                           weight, a_prev,
                                                           a_star)
        self.assertAlmostEqualVec(cape_consumed, zmtc_cape_consumed)


if __name__ == "__main__":
    unittest.main()
