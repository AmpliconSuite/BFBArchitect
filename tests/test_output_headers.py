import tempfile
import unittest
from pathlib import Path

from bfbarchitect.BFBArchitect import (
    BFB_OUTPUT_HEADER,
    write_bfb_cycles,
    write_bfb_graph,
)
from bfbarchitect.BFBVisualizer import parse_scores
from bfbarchitect.graph_input import parse_graph_file


class OutputHeaderTests(unittest.TestCase):
    def test_graph_writer_emits_version_header_and_remains_parseable(self):
        segments = [
            ("chr1", 100, 199, 4, 12.5, 7),
            ("chr1", 200, 299, 4, 12.0, 9),
        ]

        with tempfile.TemporaryDirectory() as tmpdir:
            graph_path = Path(tmpdir) / "sample_graph.txt"

            write_bfb_graph(graph_path, segments, [], {})

            lines = graph_path.read_text().splitlines()
            self.assertEqual(lines[0], BFB_OUTPUT_HEADER.rstrip())
            self.assertTrue(lines[1].startswith("SequenceEdge:"))

            svs, chrom_segs = parse_graph_file(graph_path)
            self.assertEqual(svs, [])
            self.assertEqual(len(chrom_segs["chr1"]), 2)

    def test_cycles_writer_emits_version_header_and_remains_parseable(self):
        segments = [("chr1", 100, 199, 4, 12.5, 7)]

        with tempfile.TemporaryDirectory() as tmpdir:
            cycles_path = Path(tmpdir) / "sample_cycles.txt"

            write_bfb_cycles(cycles_path, segments, [[1]], [0.25], 2)

            lines = cycles_path.read_text().splitlines()
            self.assertEqual(lines[0], BFB_OUTPUT_HEADER.rstrip())
            self.assertTrue(lines[1].startswith("Interval\t"))

            scores, seg_num = parse_scores(cycles_path)
            self.assertEqual(seg_num, 1)
            self.assertEqual(scores[0]["Final_score"], "0.25")
            self.assertEqual(scores[0]["Multiplicity"], "2")


if __name__ == "__main__":
    unittest.main()
