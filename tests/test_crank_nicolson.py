import unittest

import numpy as np

from diffusion import analytical_profile, solve_slab_cn, solve_tridiagonal


class CrankNicolsonTests(unittest.TestCase):
    def test_tridiagonal_against_dense_solver(self):
        lower = np.array([-0.5, -0.7, -0.2])
        diagonal = np.array([2.0, 2.5, 3.0, 1.5])
        upper = np.array([-0.3, -0.4, -0.6])
        rhs = np.array([1.0, 0.2, -0.3, 0.7])
        matrix = np.diag(diagonal)+np.diag(lower, -1)+np.diag(upper, 1)
        originals = [values.copy() for values in (lower, diagonal, upper, rhs)]
        np.testing.assert_allclose(solve_tridiagonal(lower, diagonal, upper, rhs), np.linalg.solve(matrix, rhs), atol=1e-14)
        for original, current in zip(originals, (lower, diagonal, upper, rhs)):
            np.testing.assert_array_equal(original, current)

    def test_matches_series(self):
        for state in solve_slab_cn([0.01, 0.05, 0.2], nodes=101, r=1):
            error = np.max(np.abs(state.concentration-analytical_profile(state.position, state.tau)))
            self.assertLess(error, 5e-4)
            np.testing.assert_allclose(state.concentration, state.concentration[::-1], atol=1e-13)
            np.testing.assert_array_equal(state.concentration[[0, -1]], [1, 1])

    def test_grid_convergence(self):
        errors = []
        for nodes in [21, 41, 81]:
            state = solve_slab_cn([0.05], nodes=nodes, r=1)[0]
            errors.append(np.max(np.abs(state.concentration-analytical_profile(state.position, state.tau))))
        self.assertLess(errors[1], errors[0]*0.35)
        self.assertLess(errors[2], errors[1]*0.35)

    def test_stability_does_not_guarantee_bounds(self):
        # One deliberately coarse step can overshoot despite linear stability.
        state = solve_slab_cn([2.5], nodes=3, r=10)[0]
        self.assertAlmostEqual(state.concentration[1], 20/11)
        self.assertGreater(state.concentration[1], 1)

    def test_times_and_initial_state(self):
        states = solve_slab_cn([0, 0.00317, 0.02741], nodes=11)
        self.assertEqual([s.tau for s in states], [0, 0.00317, 0.02741])
        np.testing.assert_array_equal(states[0].concentration, [1]+[0]*9+[1])

    def test_invalid_inputs(self):
        for times in [[], [-1], [0.1, 0.1], [np.inf]]:
            with self.assertRaises(ValueError):
                solve_slab_cn(times)
        for r in [0, -1, np.nan]:
            with self.assertRaises(ValueError):
                solve_slab_cn([0.1], r=r)
        with self.assertRaises(ValueError):
            solve_tridiagonal([], [0], [], [1])


if __name__ == '__main__':
    unittest.main()
