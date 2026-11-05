import subprocess
import sys
import unittest

from bfbarchitect.BFBSolver import reconstruct_BFB_cbc


class CbcSolverTests(unittest.TestCase):
    def test_cbc_finds_known_objectives(self):
        sequence, objective = reconstruct_BFB_cbc([2], [0], [0], 1, max_time=60, max_threads=1)
        self.assertEqual(sequence, [1, -1])
        self.assertAlmostEqual(objective, 2.0, places=6)

        sequence, objective = reconstruct_BFB_cbc([6, 4, 2], [1, 0, 1], [1, 1, 0], 1, max_time=60, max_threads=1)
        self.assertEqual(len(sequence), 12)
        self.assertAlmostEqual(objective, 4.0, places=6)

    def test_package_imports_without_pulp(self):
        code = (
            "import sys\n"
            "sys.modules['pulp'] = None\n"
            "import bfbarchitect\n"
            "from bfbarchitect.BFBSolver import reconstruct_BFB_cbc\n"
            "try:\n"
            "    reconstruct_BFB_cbc([2], [0], [0], 1)\n"
            "except ImportError as error:\n"
            "    print('pulp' in str(error))\n"
        )
        result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "True")


if __name__ == "__main__":
    unittest.main()
