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

# Constants used in every test.
tau = 2700.
deltat = 1800.
# We pick arbitrary nondimensional values for CAPE.
a_crit = 1.
# Number of cases in each category examined below.
num_cases = 3

# Inactive cases
a_prev_inactive = np.array([0.25, 1., 0.5])
a_star_inactive = np.array([1., 0.5, 0.25])

# Ongoing cases
a_prev_ongoing = np.array([1.25, 2., 1.5])
a_star_ongoing = np.array([1.5, 2., 1.4])

# Starting cases
a_prev_starting = np.array([0.25, 1., 0.5])
a_star_starting = np.array([1.3, 1.5, 2.])

# Ending cases
a_prev_ending = np.array([1.3, 1.5, 2.])
a_star_ending = np.array([0.25, 1., 1.01])

class TestCapeConsumptions(VectorCase):

    def test_ongoing_cape_consumption(self):
        """Ongoing CAPE consumption function uses the correct formula."""
        a_prev = a_prev_ongoing
        a_star = a_star_ongoing
        weight = zmtc.weight(deltat/tau)
        cape_consumed = (a_prev - a_crit) * weight / tau + \
                        (a_star - a_prev) * (1 - weight) / deltat
        zmtc_cape_consumed = zmtc.cape_consumption_ongoing(tau, a_crit, deltat,
                                                           weight, a_prev,
                                                           a_star)
        self.assertAlmostEqualVec(cape_consumed, zmtc_cape_consumed,
                                  places=DCAPE_PLACES)

    def test_starting_cape_consumption(self):
        """Starting CAPE consumption function uses the correct formula."""
        a_prev = a_prev_starting
        a_star = a_star_starting
        f = (a_crit - a_prev) / (a_star - a_prev)
        weight = zmtc.weight_1d((1-f)*deltat/tau)
        cape_consumed = (a_star - a_crit) * (1. - weight) / deltat
        zmtc_cape_consumed = zmtc.cape_consumption_starting(tau, a_crit, deltat,
                                                            a_prev, a_star)
        self.assertAlmostEqualVec(cape_consumed, zmtc_cape_consumed,
                                  places=DCAPE_PLACES)

    def test_ending_cape_consumption(self):
        """Ending CAPE consumption function uses the correct formula."""
        a_prev = a_prev_ending
        a_star = a_star_ending
        f = zmtc.end_time_frac(tau, a_crit, deltat, a_prev, a_star)
        assert all(f < 1.), "incorrect input for ending CAPE consumption"
        cape_consumed = (f * a_star + (1. - f) * a_prev - a_crit) / deltat
        zmtc_cape_consumed = zmtc.cape_consumption_ending(a_crit, deltat, f,
                                                          a_prev, a_star)
        self.assertAlmostEqualVec(cape_consumed, zmtc_cape_consumed,
                                  places=DCAPE_PLACES)


class TestCapeConsumptionRate(VectorCase):

    def test_cape_consumption_rate_inactive(self):
        """Test the combined consumption rate function on inactive cases."""
        a_prev = a_prev_inactive
        a_star = a_star_inactive
        cape_consumed = np.zeros(a_prev.shape)
        zmtc_cape_consumed = zmtc.cape_consumption_rate(tau, a_crit, deltat,
                                                        a_prev, a_star)
        self.assertAlmostEqualVec(cape_consumed, zmtc_cape_consumed,
                                  places=DCAPE_PLACES)

    def test_cape_consumption_rate_ongoing(self):
        """Test the combined consumption rate function on ongoing cases."""
        a_prev = a_prev_ongoing
        a_star = a_star_ongoing
        weight = zmtc.weight(deltat/tau)
        cape_consumed = zmtc.cape_consumption_ongoing(tau, a_crit, deltat,
                                                           weight, a_prev,
                                                           a_star)
        zmtc_cape_consumed = zmtc.cape_consumption_rate(tau, a_crit, deltat,
                                                        a_prev, a_star)
        self.assertAlmostEqualVec(cape_consumed, zmtc_cape_consumed,
                                  places=DCAPE_PLACES)

    def test_cape_consumption_rate_starting(self):
        """Test the combined consumption rate function on starting cases."""
        a_prev = a_prev_starting
        a_star = a_star_starting
        cape_consumed = zmtc.cape_consumption_starting(tau, a_crit, deltat,
                                                       a_prev, a_star)
        zmtc_cape_consumed = zmtc.cape_consumption_rate(tau, a_crit, deltat,
                                                        a_prev, a_star)
        self.assertAlmostEqualVec(cape_consumed, zmtc_cape_consumed,
                                  places=DCAPE_PLACES)

    def test_cape_consumption_rate_ending(self):
        """Test the combined consumption rate function on ending cases."""
        a_prev = a_prev_ending
        a_star = a_star_ending
        f = zmtc.end_time_frac(tau, a_crit, deltat, a_prev, a_star)
        cape_consumed = zmtc.cape_consumption_ending(a_crit, deltat, f,
                                                     a_prev, a_star)
        zmtc_cape_consumed = zmtc.cape_consumption_rate(tau, a_crit, deltat,
                                                        a_prev, a_star)
        self.assertAlmostEqualVec(cape_consumed, zmtc_cape_consumed,
                                  places=DCAPE_PLACES)

    def test_cape_consumption_rate_mixed(self):
        """The combined consumption rate correctly identifies mixed cases."""
        # These arrays will contain all the cases above jumbled up.
        a_prev = np.zeros((4*num_cases,))
        a_star = np.zeros((4*num_cases,))
        cape_consumed = np.zeros((4*num_cases,))
        # Inactive cases.
        cape_inactive = zmtc.cape_consumption_rate(tau, a_crit, deltat,
                                                   a_prev_inactive,
                                                   a_star_inactive)
        for i in range(num_cases):
            big_i = 4*i
            a_prev[big_i] = a_prev_inactive[i]
            a_star[big_i] = a_star_inactive[i]
            cape_consumed[big_i] = cape_inactive[i]
        # Ongoing cases.
        cape_ongoing = zmtc.cape_consumption_rate(tau, a_crit, deltat,
                                                  a_prev_ongoing,
                                                  a_star_ongoing)
        for i in range(num_cases):
            big_i = 4*i+1
            a_prev[big_i] = a_prev_ongoing[i]
            a_star[big_i] = a_star_ongoing[i]
            cape_consumed[big_i] = cape_ongoing[i]
        # Starting cases.
        cape_starting = zmtc.cape_consumption_rate(tau, a_crit, deltat,
                                                   a_prev_starting,
                                                   a_star_starting)
        for i in range(num_cases):
            big_i = 4*i+2
            a_prev[big_i] = a_prev_starting[i]
            a_star[big_i] = a_star_starting[i]
            cape_consumed[big_i] = cape_starting[i]
        # Ending cases.
        cape_ending = zmtc.cape_consumption_rate(tau, a_crit, deltat,
                                                   a_prev_ending,
                                                   a_star_ending)
        for i in range(num_cases):
            big_i = 4*i+3
            a_prev[big_i] = a_prev_ending[i]
            a_star[big_i] = a_star_ending[i]
            cape_consumed[big_i] = cape_ending[i]
        zmtc_cape_consumed = zmtc.cape_consumption_rate(tau, a_crit, deltat,
                                                        a_prev, a_star)
        self.assertAlmostEqualVec(cape_consumed, zmtc_cape_consumed,
                                  places=DCAPE_PLACES)



if __name__ == "__main__":
    unittest.main()
