import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1] / "bfbarchitect"


class ImportModeTests(unittest.TestCase):
    def test_modules_run_directly_as_scripts(self):
        for script in ("BFBArchitect.py", "BFBVisualizer.py", "SVCaller.py", "batch_run.py"):
            with self.subTest(script=script), tempfile.TemporaryDirectory() as tmpdir:
                result = subprocess.run(
                    [sys.executable, str(PACKAGE_DIR / script), "--help"],
                    cwd=tmpdir, capture_output=True, text=True,
                )
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_package_import_reports_the_real_error(self):
        code = "import sys\nsys.modules['pysam'] = None\nimport bfbarchitect\n"
        result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("pysam", result.stderr)
        self.assertNotIn("_version", result.stderr)

    def test_direct_run_reports_the_real_error(self):
        code = (
            "import runpy, sys\n"
            "sys.modules['pysam'] = None\n"
            "sys.path.insert(0, sys.argv[1])\n"
            "sys.argv = ['BFBArchitect.py', '--help']\n"
            "runpy.run_path(sys.path[0] + '/BFBArchitect.py', run_name='__main__')\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            result = subprocess.run(
                [sys.executable, "-c", code, str(PACKAGE_DIR)],
                cwd=tmpdir, capture_output=True, text=True,
            )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("pysam", result.stderr)


if __name__ == "__main__":
    unittest.main()
