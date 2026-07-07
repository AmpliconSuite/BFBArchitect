import logging
import unittest
from unittest.mock import patch

import bfbarchitect.BFBArchitect as bfb


class SolverSelectionTests(unittest.TestCase):
    def test_detect_solver_validates_gurobi_before_selecting_it(self):
        with patch.object(bfb, "_gurobi_availability_error", return_value=None), \
             patch.object(bfb, "_mosek_available", return_value=False):
            self.assertEqual(bfb.detect_solver(), "gurobi")

    def test_detect_solver_falls_back_when_gurobi_license_checkout_fails(self):
        with patch.object(bfb, "_gurobi_availability_error",
                          return_value="Could not resolve host: token.gurobi.com"), \
             patch.object(bfb, "_mosek_available", return_value=False):
            self.assertEqual(bfb.detect_solver(), "cbc")

    def test_multiple_with_non_gurobi_warns_and_runs_single_solution(self):
        with patch.object(bfb, "reconstruct_BFB_cbc", return_value=([1], 0.0)) as cbc:
            with self.assertLogs("BFBArchitect", level=logging.WARNING) as logs:
                strings, obj_val = bfb._solve_bfb_ilp(
                    [1], [0], [0], 1,
                    solver="cbc", multiple=True, auto_solver=False)

        self.assertEqual(strings, [[1]])
        self.assertEqual(obj_val, 0.0)
        cbc.assert_called_once()
        self.assertIn("--multiple requires Gurobi", logs.output[0])

    def test_auto_selected_gurobi_retries_cbc_on_token_failure(self):
        with patch.object(
            bfb,
            "reconstruct_BFB_gurobi",
            side_effect=RuntimeError("Could not resolve host: token.gurobi.com"),
        ) as gurobi, \
             patch.object(bfb, "_detect_non_gurobi_solver", return_value="cbc"), \
             patch.object(bfb, "reconstruct_BFB_cbc", return_value=([1], 0.0)) as cbc:
            with self.assertLogs("BFBArchitect", level=logging.WARNING) as logs:
                strings, obj_val = bfb._solve_bfb_ilp(
                    [1], [0], [0], 1,
                    solver="gurobi", multiple=False, auto_solver=True)

        self.assertEqual(strings, [[1]])
        self.assertEqual(obj_val, 0.0)
        gurobi.assert_called_once()
        cbc.assert_called_once()
        self.assertIn("retrying with CBC solver", logs.output[0])

    def test_explicit_gurobi_does_not_fallback_on_token_failure(self):
        with patch.object(
            bfb,
            "reconstruct_BFB_gurobi",
            side_effect=RuntimeError("Could not resolve host: token.gurobi.com"),
        ), \
             patch.object(bfb, "reconstruct_BFB_cbc") as cbc:
            with self.assertRaisesRegex(RuntimeError, "token.gurobi.com"):
                bfb._solve_bfb_ilp(
                    [1], [0], [0], 1,
                    solver="gurobi", multiple=False, auto_solver=False)

        cbc.assert_not_called()


if __name__ == "__main__":
    unittest.main()
