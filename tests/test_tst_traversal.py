import unittest
from unittest.mock import patch

from bfbarchitect.datatypes import SV
import bfbarchitect.graph_input as graph_input
from bfbarchitect.graph_input import (
    _build_graph_lookups,
    _dedupe_reciprocal_local_rescue_anchors,
    _directional_tst_pair_candidates,
    _filter_illegal_nested_foldbacks,
    whole_graph_as_region,
)


class DirectionalTstTraversalTests(unittest.TestCase):
    def _entries_and_segments(self, sv2_shard_strand="+"):
        sv1 = SV("chr1", 5000, "-", "chr2", 2000, "+")
        sv2 = SV("chr1", 5300, sv2_shard_strand, "chr2", 2101, "+")
        entries = [(sv1, 3.0, 10), (sv2, 3.0, 10)]
        chrom_segs = {
            "chr1": [
                (5000, 5100, 4.0, 0.0, 0),
                (5101, 5300, 4.0, 0.0, 0),
            ],
            "chr2": [
                (1000, 2000, 10.0, 0.0, 0),
                (2101, 3000, 5.0, 0.0, 0),
            ],
        }
        return entries, chrom_segs

    def test_two_hop_shard_path_reaches_expected_endpoint(self):
        entries, chrom_segs = self._entries_and_segments()
        bp_to_seg = _build_graph_lookups(entries, chrom_segs)

        candidates = _directional_tst_pair_candidates(
            entries[0], entries[1], bp_to_seg, chrom_segs,
            fb_dist=50000, far_min=100, shard_max_bp=5000, max_hops=5,
            region_by_chrom={"chr2": [(1, 10000)]},
        )

        chr2_candidates = [
            candidate for candidate in candidates
            if candidate["synth"].chrom1 == "chr2"
        ]
        self.assertEqual(len(chr2_candidates), 1)
        candidate = chr2_candidates[0]
        self.assertEqual(str(candidate["synth"]), "chr2:2000+->chr2:2101+")
        self.assertEqual(
            candidate["path"],
            [("chr1", 5000, 5100), ("chr1", 5101, 5300)],
        )

    def test_wrong_shard_endpoint_orientation_does_not_close_path(self):
        entries, chrom_segs = self._entries_and_segments(sv2_shard_strand="-")
        bp_to_seg = _build_graph_lookups(entries, chrom_segs)

        candidates = _directional_tst_pair_candidates(
            entries[0], entries[1], bp_to_seg, chrom_segs,
            fb_dist=50000, far_min=100, shard_max_bp=5000, max_hops=5,
            region_by_chrom={"chr2": [(1, 10000)]},
        )

        self.assertFalse([
            candidate for candidate in candidates
            if candidate["synth"].chrom1 == "chr2"
        ])

    def test_reciprocal_local_rescue_prefers_larger_core_mass(self):
        bridge_sv = object()
        weak_foldback = object()
        strong_foldback = object()
        chrom_segs = {
            "chr1": [
                (62466092, 62478528, 11.149, 0.0, 0),
                (62478529, 62491907, 12.586, 0.0, 0),
                (62491908, 62505286, 16.951, 0.0, 0),
                (62505287, 70935033, 1.355, 0.0, 0),
                (70935034, 70976824, 1.323, 0.0, 0),
                (70976825, 70977737, 7.003, 0.0, 0),
                (70977738, 71006962, 17.177, 0.0, 0),
                (71006963, 71108889, 1.347, 0.0, 0),
            ],
        }
        anchors = [
            {
                "boundary_side": "right",
                "boundary_endpoint": ("chr1", 62505286, "+"),
                "landed_endpoint": ("chr1", 71006962, "+"),
                "boundary_sv": bridge_sv,
                "foldback": strong_foldback,
                "count": 3,
                "boundary_cn": 17.198,
                "boundary_outside_cn": 1.368,
                "landed_cn": 17.177,
                "landed_neighbor_cn": 7.003,
            },
            {
                "boundary_side": "right",
                "boundary_endpoint": ("chr1", 71006962, "+"),
                "landed_endpoint": ("chr1", 62505286, "+"),
                "boundary_sv": bridge_sv,
                "foldback": weak_foldback,
                "count": 6,
                "boundary_cn": 17.177,
                "boundary_outside_cn": 1.347,
                "landed_cn": 17.198,
                "landed_neighbor_cn": 14.095,
            },
        ]

        kept = _dedupe_reciprocal_local_rescue_anchors(anchors, chrom_segs)

        self.assertEqual(len(kept), 1)
        self.assertIs(kept[0]["foldback"], strong_foldback)
        self.assertEqual(kept[0]["boundary_endpoint"], ("chr1", 62505286, "+"))

    def test_whole_graph_uses_largest_chrom_with_native_foldback(self):
        chr2_foldback = SV("chr2", 20162013, "+", "chr2", 20162804, "+")
        chr8_jump = SV("chr3", 50043258, "-", "chr8", 127747325, "+")
        svs = [
            (chr2_foldback, 2.7, 65),
            (chr8_jump, 1.2, 63),
        ]
        chrom_segs = {
            "chr2": [
                (19910152, 20008469, 4.0, 0.0, 0),
                (20008470, 20058152, 5.9, 0.0, 0),
                (20058153, 20162013, 7.1, 0.0, 0),
                (20162014, 20162804, 4.5, 0.0, 0),
                (20162805, 20265154, 1.8, 0.0, 0),
            ],
            "chr8": [
                (127188812, 127305705, 5.7, 0.0, 0),
                (127305706, 127467019, 6.8, 0.0, 0),
                (127467020, 127747325, 5.6, 0.0, 0),
            ],
        }
        region_data = ([("chr2", 19910152, 20265154, 6.0, 0.0, 0)],
                       [5], [0], [3], [chr2_foldback], {chr2_foldback: (2.7, 65)})

        with patch.object(graph_input, "parse_graph_file", return_value=(svs, chrom_segs)), \
                patch.object(graph_input, "subsect_graph_for_region",
                             return_value=[region_data]) as subsect:
            result = whole_graph_as_region("mock_graph.txt", deletion=False)

        self.assertEqual(result[-1], "chr2")
        subsect.assert_called_once()
        self.assertEqual(subsect.call_args.args[1], [("chr2", 19910152, 20265154)])

    def test_nested_opposite_polarity_foldbacks_are_illegal(self):
        outer = SV("chr12", 71315480, "-", "chr12", 71316538, "-")
        inner = SV("chr12", 71315763, "+", "chr12", 71316536, "+")

        self.assertEqual(_filter_illegal_nested_foldbacks([outer, inner]), [])

    def test_non_nested_foldback_pairs_are_kept(self):
        left = SV("chr12", 1000, "-", "chr12", 1500, "-")
        right = SV("chr12", 2000, "+", "chr12", 2500, "+")
        same_polarity_nested = SV("chr12", 1100, "-", "chr12", 1400, "-")

        kept = _filter_illegal_nested_foldbacks(
            [left, right, same_polarity_nested]
        )

        self.assertEqual(kept, [left, right, same_polarity_nested])


if __name__ == "__main__":
    unittest.main()
