"""Checks against the analytical model and physical bounds, not measured data."""

import unittest

import numpy as np

from diffusion import (
    analytical_profile, analytical_uptake, physical_time, solve_slab,
    time_for_uptake,
)


class DiffusionTests(unittest.TestCase):
    def test_matches_series(self):
        # Resolve the early boundary layer with 201 nodes.
        for state in solve_slab([0.01, 0.05, 0.2], nodes=201):
            error = np.max(np.abs(state.concentration - analytical_profile(state.position, state.tau)))
            self.assertLess(error, 5e-4)
            self.assertAlmostEqual(state.uptake, analytical_uptake(state.tau), delta=5e-4)

    def test_bounds_symmetry_and_monotone_uptake(self):
        states = solve_slab([0, 0.002, 0.01, 0.1, 1], nodes=51)
        for state in states:
            self.assertGreaterEqual(state.concentration.min(), 0)
            self.assertLessEqual(state.concentration.max(), 1)
            np.testing.assert_allclose(state.concentration, state.concentration[::-1], atol=1e-14)
            np.testing.assert_array_equal(state.concentration[[0, -1]], [1, 1])
        for a, b in zip(states, states[1:]):
            self.assertGreater(b.uptake, a.uptake)
            self.assertTrue(np.all(b.concentration >= a.concentration - 1e-14))
        self.assertGreater(states[-1].concentration.min(), 0.9999)

    def test_grid_refinement(self):
        errors = []
        for nodes in [21, 41, 81, 161]:
            state = solve_slab([0.05], nodes=nodes)[0]
            errors.append(np.max(np.abs(state.concentration - analytical_profile(state.position, state.tau))))
        for coarse, fine in zip(errors, errors[1:]):
            self.assertLess(fine, 0.35 * coarse)

    def test_requested_output_times(self):
        times = [0.0, 0.00317, 0.02741]
        self.assertEqual([state.tau for state in solve_slab(times)], times)

    def test_physical_scaling(self):
        base = physical_time(0.2, 1e-6, 1e-14)
        self.assertAlmostEqual(float(base), 20)
        self.assertAlmostEqual(float(physical_time(0.2, 2e-6, 1e-14)), 4 * base)
        self.assertAlmostEqual(float(physical_time(0.2, 1e-6, 2e-14)), base / 2)

    def test_series_initial_and_long_time_limits(self):
        z = np.linspace(0, 1, 11)
        np.testing.assert_array_equal(analytical_profile(z, 0), [1]+[0]*9+[1])
        np.testing.assert_allclose(analytical_profile(z, 5), np.ones(11))
        self.assertEqual(analytical_uptake(0), 0)
        self.assertAlmostEqual(analytical_uptake(5), 1)

    def test_ninety_percent_uptake(self):
        tau = time_for_uptake(0.9)
        self.assertAlmostEqual(analytical_uptake(tau), 0.9, places=12)
        self.assertGreater(tau, 0.2)
        self.assertLess(tau, 0.22)

    def test_invalid_solver_inputs(self):
        for times in [[], [-1], [0.2, 0.1], [0.1, 0.1], [np.nan], [np.inf], [[0.1]]]:
            with self.subTest(times=times), self.assertRaises(ValueError):
                solve_slab(times)
        for nodes in [2, 3.5, True]:
            with self.subTest(nodes=nodes), self.assertRaises(ValueError):
                solve_slab([0.1], nodes=nodes)
        for r in [0, -1, 0.51, np.inf, np.nan]:
            with self.subTest(r=r), self.assertRaises(ValueError):
                solve_slab([0.1], r=r)

    def test_invalid_physical_inputs(self):
        for thickness, diffusivity in [(0, 1e-14), (1e-6, -1), (np.nan, 1)]:
            with self.assertRaises(ValueError):
                physical_time(0.1, thickness, diffusivity)

    def test_unsupported_uptake_fraction(self):
        for fraction in [0, 0.001, 1, np.nan]:
            with self.assertRaises(ValueError):
                time_for_uptake(fraction)


if __name__ == '__main__':
    unittest.main()
