import contextlib
import io
import logging
import sys
import types
import unittest
from unittest.mock import patch

import pulp
import bfbarchitect.BFBSolver as solver
import bfbarchitect.BFBArchitect as bfb
import bfbarchitect.batch_run as batch
from bfbarchitect.BFBArchitect import compute_bfb_scores


class CbcTieBreakTests(unittest.TestCase):
    cn, lf, rf = [2, 3], [0, 0], [0, 1]

    def test_default_canonicalizes_tied_solutions(self):
        with patch.object(solver, '_canonical_cbc_sequences',
                          wraps=solver._canonical_cbc_sequences) as canonical:
            sequence, objective = solver.reconstruct_BFB_cbc(
                self.cn, self.lf, self.rf, 1, max_time=30, max_threads=1)
        canonical.assert_called_once()
        self.assertEqual(sequence, [1, 2, -2, -1, 1, 2])
        self.assertAlmostEqual(objective, 3.0)

    def test_public_reconstruction_is_deterministic_by_default(self):
        segments = [('chr1', 1, 100, 3, 0, 0), ('chr1', 101, 200, 4, 0, 0)]
        strings, scores, multiplicity = bfb.reconstruct_bfb(
            segments, self.cn, self.lf, self.rf, 0, solver='cbc', silent=True, threads=1)
        self.assertEqual(strings, [[1, 2, -2, -1, 1, 2]])
        self.assertAlmostEqual(scores[0], 10 / 3)
        self.assertEqual(multiplicity, 1)

    def test_cli_has_no_tiebreak_switch(self):
        output = io.StringIO()
        with patch.object(sys, 'argv', ['BFBArchitect.py', '--help']), \
                contextlib.redirect_stdout(output), self.assertRaises(SystemExit) as result:
            bfb.main()
        self.assertEqual(result.exception.code, 0)
        self.assertNotIn('--deterministic-cbc', output.getvalue())

    def test_batch_cli_has_no_tiebreak_switch(self):
        self.assertNotIn('--deterministic-cbc', batch.build_parser().format_help())

    def test_automatic_fallback_to_cbc_is_deterministic(self):
        with patch.object(bfb, 'reconstruct_BFB_gurobi',
                          side_effect=RuntimeError('Gurobi license unavailable')), \
                patch.object(bfb, '_detect_non_gurobi_solver', return_value='cbc'):
            strings, objective = bfb._solve_bfb_ilp(
                self.cn, self.lf, self.rf, 1, solver='gurobi', multiple=False,
                auto_solver=True, warn=False, threads=1)
        self.assertEqual(strings, [[1, 2, -2, -1, 1, 2]])
        self.assertAlmostEqual(objective, 3.0)

    def test_tied_solutions_are_identical_across_seeds_and_threads(self):
        command = solver._cbc_command
        expected = {}
        for start in (1, -1):
            for seed, threads in ((1, 1), (42, 1), (7, 2)):
                def seeded_command(max_time, max_threads):
                    result = command(max_time, max_threads)
                    result.options = [f'RandomS {seed}', f'randomCbcSeed {seed}']
                    return result
                with self.subTest(start=start, seed=seed, threads=threads), \
                        patch.object(solver, '_cbc_command', side_effect=seeded_command):
                    sequence, objective = solver.reconstruct_BFB_cbc(
                        self.cn, self.lf, self.rf, start, max_time=30, max_threads=threads)
                    self.assertAlmostEqual(objective, 3.0 if start == 1 else 2.0)
                    self.assertTrue(solver.check_BFB_string(sequence))
                    score = compute_bfb_scores(self.cn, self.lf, self.rf, [sequence], 1,
                                               logging.getLogger('test'), silent=True)
                    if start not in expected:
                        expected[start] = (sequence, score)
                    self.assertEqual((sequence, score), expected[start])
        self.assertEqual(expected[1][0], [1, 2, -2, -1, 1, 2])
        self.assertAlmostEqual(expected[1][1][0], 10 / 3)

    def test_primary_objective_is_preserved(self):
        original_solve = pulp.LpProblem.solve
        for start in (1, -1):
            primary = []

            def record_primary(model, command):
                result = original_solve(model, command)
                if not primary:
                    primary.append(model.objective.value())
                return result

            with self.subTest(start=start), patch.object(pulp.LpProblem, 'solve', new=record_primary):
                _, objective = solver.reconstruct_BFB_cbc(
                    self.cn, self.lf, self.rf, start, max_time=30, max_threads=1)
                self.assertAlmostEqual(objective, primary[0])
                self.assertAlmostEqual(objective, 3.0 if start == 1 else 2.0)

    def test_interrupted_tiebreak_keeps_the_original_feasible_sequence(self):
        original_solve = pulp.LpProblem.solve
        snapshots = []

        def interrupted(model, command):
            if snapshots:
                return types.SimpleNamespace(has_solution=False,
                                             status=types.SimpleNamespace(name='TimeLimit'))
            result = original_solve(model, command)
            sequences = []
            for var in model.variables():
                if var.name.startswith('cs_') and var.value() > 0.5:
                    sequences.append(tuple(map(int, var.name[3:].split('_'))))
            snapshots.append(solver._build_bfb_string(sequences))
            return result

        with patch.object(pulp.LpProblem, 'solve', new=interrupted), \
                self.assertLogs('BFBArchitect', level='WARNING'):
            sequence, objective = solver.reconstruct_BFB_cbc(
                self.cn, self.lf, self.rf, 1, max_time=30, max_threads=1)
        self.assertEqual(sequence, snapshots[0])
        self.assertAlmostEqual(objective, 3.0)

    def test_time_limited_primary_solution_skips_tiebreak(self):
        original_solve = pulp.LpProblem.solve

        def time_limited(model, command):
            original_solve(model, command)
            return types.SimpleNamespace(has_solution=True,
                                         status=types.SimpleNamespace(name='TimeLimit'))

        with patch.object(pulp.LpProblem, 'solve', new=time_limited), \
                patch.object(solver, '_canonical_cbc_sequences') as canonical, \
                self.assertLogs('BFBArchitect', level='WARNING'):
            sequence, _ = solver.reconstruct_BFB_cbc(
                self.cn, self.lf, self.rf, 1, max_time=30, max_threads=1)
        canonical.assert_not_called()
        self.assertTrue(solver.check_BFB_string(sequence))

    def test_exhausted_budget_keeps_incumbent_without_more_solves(self):
        original_solve = pulp.LpProblem.solve
        snapshots = []

        def solve_once(model, command):
            result = original_solve(model, command)
            snapshots.append(result)
            return result

        with patch.object(pulp.LpProblem, 'solve', new=solve_once), \
                patch.object(solver._time, 'monotonic', side_effect=[0, 31]), \
                self.assertLogs('BFBArchitect', level='WARNING'):
            sequence, objective = solver.reconstruct_BFB_cbc(
                self.cn, self.lf, self.rf, 1, max_time=30, max_threads=1)
        self.assertEqual(len(snapshots), 1)
        self.assertTrue(solver.check_BFB_string(sequence))
        self.assertAlmostEqual(objective, 3.0)

    def test_no_integer_solution_raises_instead_of_returning_a_fractional_string(self):
        result = types.SimpleNamespace(has_solution=False)
        with patch.object(pulp.LpProblem, 'solve', return_value=result):
            with self.assertRaisesRegex(RuntimeError, 'feasible integer BFB solution'):
                solver.reconstruct_BFB_cbc(self.cn, self.lf, self.rf, 1)

    def test_legacy_time_limit_is_not_treated_as_proven_optimal(self):
        model = types.SimpleNamespace(sol_status=2)
        self.assertTrue(solver._cbc_has_solution(model, 1))
        self.assertFalse(solver._cbc_is_optimal(model, 1))


if __name__ == '__main__':
    unittest.main()
