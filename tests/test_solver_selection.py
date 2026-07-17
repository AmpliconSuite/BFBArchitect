import logging
import os
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

import bfbarchitect.BFBArchitect as bfb


class SolverSelectionTests(unittest.TestCase):
    @staticmethod
    def _fake_gurobi_module(expected_license=None, home_default=False):
        module = types.ModuleType("gurobipy")

        class Env:
            def __init__(self, empty=False):
                self.empty = empty

            def setParam(self, name, value):
                pass

            def start(self):
                if expected_license is not None:
                    if os.environ.get("GRB_LICENSE_FILE") != expected_license:
                        raise RuntimeError("No Gurobi license found")
                elif home_default:
                    if not Path(os.path.expanduser("~/gurobi.lic")).is_file():
                        raise RuntimeError("No Gurobi license found")

            def dispose(self):
                pass

        module.Env = Env
        return module

    def test_gurobi_selected_through_grb_license_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            license_file = Path(tmpdir) / "mounted" / "gurobi.lic"
            license_file.parent.mkdir()
            license_file.touch()
            fake_gurobi = self._fake_gurobi_module(str(license_file))

            with patch.dict(os.environ, {
                    "HOME": str(Path(tmpdir) / "root"),
                    "GRB_LICENSE_FILE": str(license_file),
            }, clear=True), patch.dict(sys.modules, {"gurobipy": fake_gurobi}):
                self.assertEqual(bfb.detect_solver(), "gurobi")

    def test_mosek_selected_through_directory_in_license_variable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            license_dir = Path(tmpdir) / "mounted-mosek"
            license_dir.mkdir()
            (license_dir / "mosek.lic").touch()
            fake_mosek = types.ModuleType("mosek")

            with patch.dict(os.environ, {
                    "HOME": str(Path(tmpdir) / "root"),
                    "MOSEKLM_LICENSE_FILE": str(license_dir),
            }, clear=True), patch.dict(sys.modules, {"mosek": fake_mosek}), \
                    patch.object(bfb, "_gurobi_availability_error",
                                 return_value="gurobipy is not installed"):
                self.assertEqual(bfb.detect_solver(), "mosek")

    def test_mosek_license_variable_accepts_documented_formats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            license_file = Path(tmpdir) / "mosek.lic"
            license_file.touch()
            values = (
                str(license_file),
                "@licenses.example.org",
                "27007@licenses.example.org",
                f"{Path(tmpdir) / 'missing.lic'}{os.pathsep}{license_file}",
                "START_LICENSE\nFEATURE PTS redacted\nEND_LICENSE",
            )
            for value in values:
                with self.subTest(value=value.splitlines()[0]), \
                        patch.dict(os.environ, {"MOSEKLM_LICENSE_FILE": value}, clear=True):
                    self.assertTrue(bfb._mosek_license_configured())

    def test_home_license_locations_remain_supported(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            home = Path(tmpdir)
            (home / "gurobi.lic").touch()
            (home / "mosek").mkdir()
            (home / "mosek" / "mosek.lic").touch()
            fake_gurobi = self._fake_gurobi_module(home_default=True)
            fake_mosek = types.ModuleType("mosek")

            with patch.dict(os.environ, {"HOME": str(home)}, clear=True), \
                    patch.dict(sys.modules, {
                        "gurobipy": fake_gurobi,
                        "mosek": fake_mosek,
                    }):
                self.assertIsNone(bfb._gurobi_availability_error())
                self.assertTrue(bfb._mosek_available())

    def test_detect_solver_validates_gurobi_before_selecting_it(self):
        with patch.object(bfb, "_gurobi_availability_error", return_value=None), \
             patch.object(bfb, "_mosek_available", return_value=False):
            self.assertEqual(bfb.detect_solver(), "gurobi")

    def test_detect_solver_falls_back_when_gurobi_license_checkout_fails(self):
        with patch.object(bfb, "_gurobi_availability_error",
                          return_value="Could not resolve host: token.gurobi.com"), \
             patch.object(bfb, "_mosek_available", return_value=False):
            self.assertEqual(bfb.detect_solver(), "cbc")

    def test_missing_solver_packages_select_cbc(self):
        with patch.dict(os.environ, {"HOME": "/nonexistent"}, clear=True), \
                patch.dict(sys.modules, {"gurobipy": None, "mosek": None}):
            self.assertEqual(bfb.detect_solver(), "cbc")

    def test_installed_packages_without_licenses_select_cbc(self):
        fake_gurobi = self._fake_gurobi_module(expected_license="unconfigured")
        fake_mosek = types.ModuleType("mosek")
        with tempfile.TemporaryDirectory() as tmpdir, \
                patch.dict(os.environ, {"HOME": tmpdir}, clear=True), \
                patch.dict(sys.modules, {
                    "gurobipy": fake_gurobi,
                    "mosek": fake_mosek,
                }):
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

    def test_auto_selected_gurobi_retries_mosek(self):
        with patch.object(
            bfb,
            "reconstruct_BFB_gurobi",
            side_effect=RuntimeError("WLS network is unreachable"),
        ) as gurobi, \
             patch.object(bfb, "_detect_non_gurobi_solver", return_value="mosek"), \
             patch.object(bfb, "reconstruct_BFB_mosek", return_value=([[1]], 0.0)) as mosek, \
             patch.object(bfb, "reconstruct_BFB_cbc") as cbc:
            with self.assertLogs("BFBArchitect", level=logging.WARNING) as logs:
                strings, obj_val = bfb._solve_bfb_ilp(
                    [1], [0], [0], 1,
                    solver="gurobi", multiple=False, auto_solver=True)

        self.assertEqual(strings, [[1]])
        self.assertEqual(obj_val, 0.0)
        gurobi.assert_called_once()
        mosek.assert_called_once()
        cbc.assert_not_called()
        self.assertIn("retrying with MOSEK solver", logs.output[0])

    def test_auto_selection_retries_cbc_when_mosek_license_fails(self):
        with patch.object(
            bfb,
            "reconstruct_BFB_gurobi",
            side_effect=RuntimeError("WLS network is unreachable"),
        ), patch.object(bfb, "_detect_non_gurobi_solver", return_value="mosek"), \
             patch.object(
                 bfb,
                 "reconstruct_BFB_mosek",
                 side_effect=RuntimeError("MOSEK license server is not responding"),
             ) as mosek, \
             patch.object(bfb, "reconstruct_BFB_cbc", return_value=([1], 0.0)) as cbc:
            with self.assertLogs("BFBArchitect", level=logging.WARNING) as logs:
                strings, obj_val = bfb._solve_bfb_ilp(
                    [1], [0], [0], 1,
                    solver="gurobi", multiple=False, auto_solver=True)

        self.assertEqual(strings, [[1]])
        self.assertEqual(obj_val, 0.0)
        mosek.assert_called_once()
        cbc.assert_called_once()
        self.assertEqual(len(logs.output), 2)
        self.assertIn("retrying with MOSEK solver", logs.output[0])
        self.assertIn("retrying with CBC solver", logs.output[1])

    def test_gurobi_network_error_code_triggers_automatic_fallback(self):
        error = RuntimeError("remote service unavailable")
        error.errno = 10022
        with patch.object(bfb, "reconstruct_BFB_gurobi", side_effect=error), \
             patch.object(bfb, "_detect_non_gurobi_solver", return_value="cbc"), \
             patch.object(bfb, "reconstruct_BFB_cbc", return_value=([1], 0.0)) as cbc:
            bfb._solve_bfb_ilp(
                [1], [0], [0], 1,
                solver="gurobi", multiple=False, auto_solver=True, warn=False)
        cbc.assert_called_once()

    def test_non_availability_solver_error_is_not_retried(self):
        with patch.object(
            bfb,
            "reconstruct_BFB_gurobi",
            side_effect=RuntimeError("invalid model construction"),
        ), patch.object(bfb, "reconstruct_BFB_cbc") as cbc:
            with self.assertRaisesRegex(RuntimeError, "invalid model construction"):
                bfb._solve_bfb_ilp(
                    [1], [0], [0], 1,
                    solver="gurobi", multiple=False, auto_solver=True)
        cbc.assert_not_called()

    def test_inline_license_contents_are_redacted_from_fallback_log(self):
        license_text = "START_LICENSE\nFEATURE PTS highly-secret\nEND_LICENSE"
        with patch.dict(os.environ, {"MOSEKLM_LICENSE_FILE": license_text}, clear=True), \
             patch.object(
                 bfb,
                 "reconstruct_BFB_mosek",
                 side_effect=RuntimeError(f"Rejected {license_text}"),
             ), patch.object(bfb, "reconstruct_BFB_cbc", return_value=([1], 0.0)):
            with self.assertLogs("BFBArchitect", level=logging.WARNING) as logs:
                bfb._solve_bfb_ilp(
                    [1], [0], [0], 1,
                    solver="mosek", multiple=False, auto_solver=True)

        combined_logs = "\n".join(logs.output)
        self.assertNotIn("highly-secret", combined_logs)
        self.assertIn("details redacted", combined_logs)

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

    def test_explicit_mosek_does_not_fallback_on_license_failure(self):
        with patch.object(
            bfb,
            "reconstruct_BFB_mosek",
            side_effect=RuntimeError("MOSEK license server is not responding"),
        ), patch.object(bfb, "reconstruct_BFB_cbc") as cbc:
            with self.assertRaisesRegex(RuntimeError, "license server"):
                bfb._solve_bfb_ilp(
                    [1], [0], [0], 1,
                    solver="mosek", multiple=False, auto_solver=False)

        cbc.assert_not_called()


if __name__ == "__main__":
    unittest.main()
