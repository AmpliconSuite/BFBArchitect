import itertools
import logging
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib.pyplot as plt

import bfbarchitect.BFBVisualizer as viz


class VisualizerChromosomeTests(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tmpdir.name)
        self.fonts = plt.rc_context({'font.family': 'DejaVu Sans'})
        self.fonts.__enter__()

    @staticmethod
    def close_plot_log():
        logger = logging.getLogger('BFBVisualizer')
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
            handler.close()

    def tearDown(self):
        self.close_plot_log()
        plt.close('all')
        self.fonts.__exit__(None, None, None)
        self.tmpdir.cleanup()

    def fixtures(self, graph_chrom, cnr_chrom=None, segments=1, gene_chrom=None):
        self.close_plot_log()
        start, end = 1000000, 1000000 + segments * 100000
        graph = self.root / 'sample_graph.txt'
        lines = [f'sequence\t{graph_chrom}:{start + i * 100000}-\t'
                 f'{graph_chrom}:{start + (i + 1) * 100000}+\t4\t20\t100000\t10'
                 for i in range(segments)]
        lines.append(f'discordant\t{graph_chrom}:{end - 10}+->{graph_chrom}:{end}+\t2\t10')
        graph.write_text('\n'.join(lines) + '\n')
        cycles = self.root / 'sample_cycles.txt'
        sequence = '1+,1-' if segments == 1 else '1+,2+,2-,1-'
        cycles.write_text(f'Path=1;Copy_count=1;Segments={sequence};'
                          'Path_constraints_satisfied=;Score=0;Multiplicity=1\n')
        cnr = None
        if cnr_chrom is not None:
            cnr = self.root / 'sample.cnr'
            cnr.write_text('chromosome\tstart\tend\tlog2\n'
                           f'{cnr_chrom}\t{start}\t{end}\t1\n')
        gene = None
        if gene_chrom is not None:
            gene = self.root / 'sample.gtf'
            gene.write_text(f'{gene_chrom}\tsource\tgene\t1020000\t1080000\t.\t+\t.\t'
                            'gene_id "g1"; gene_name "GENE1";\n')
        return graph, cycles, cnr, gene, start, end

    def test_foldbacks_match_both_prefix_styles(self):
        for graph_chrom, requested in itertools.product(('7', 'chr7'), repeat=2):
            with self.subTest(graph=graph_chrom, requested=requested):
                graph, _, _, _, start, end = self.fixtures(graph_chrom)
                coordinates = viz.parse_segment_coordinates(graph, 1)
                foldbacks = viz.parse_foldback_coordinate(graph, coordinates, requested, start, end)
                self.assertEqual(len(foldbacks), 1)
                self.assertEqual(foldbacks[0]['segment_name'], 'A')
                self.assertEqual(foldbacks[0]['chrom'], graph_chrom)
                self.assertEqual(coordinates['A']['chrom'], graph_chrom)

    def test_cnr_numeric_chromosomes_match_aliases_and_exclude_other_chromosomes(self):
        for cnr_chrom, requested in itertools.product(('7', 'chr7'), repeat=2):
            with self.subTest(cnr=cnr_chrom, requested=requested):
                _, _, cnr, _, start, end = self.fixtures('7', cnr_chrom)
                with cnr.open('a') as handle:
                    handle.write(f'17\t{start}\t{end}\t3\n')
                x, ranges, y = viz.extract_fcna(cnr, requested, start, end)
                self.assertEqual(x, [(start + end) // 2])
                self.assertEqual(ranges, [(start, end)])
                self.assertEqual(y, [3.0])

    def test_render_single_and_multiple_segments_with_and_without_cnr(self):
        pairs = [('7', 'chr7'), ('chr7', '7'), ('7', '7'), ('chr7', 'chr7'),
                 ('scaffold_chr7', 'scaffold_chr7')]
        for segments, (graph_chrom, annotation_chrom), with_cnr in itertools.product(
                (1, 2), pairs, (False, True)):
            with self.subTest(segments=segments, graph=graph_chrom,
                              annotation=annotation_chrom, cnr=with_cnr):
                graph, cycles, cnr, gene, start, end = self.fixtures(
                    graph_chrom, annotation_chrom if with_cnr else None,
                    segments, annotation_chrom)
                output = str(self.root / 'plot')
                with patch.object(viz, 'plot_foldbacks', wraps=viz.plot_foldbacks) as foldbacks, \
                        patch.object(viz, 'plot_structure', wraps=viz.plot_structure) as structure, \
                        patch.object(viz, 'plot_genes', wraps=viz.plot_genes) as genes:
                    viz.visualize_BFB(cycles, graph, cnr, output, gene_annotation=gene)
                self.assertTrue((self.root / 'plot_1.png').is_file())
                self.assertEqual(len(foldbacks.call_args.args[1]), 1)
                self.assertEqual(foldbacks.call_args.args[3], end - start)
                self.assertEqual(structure.call_args.args[4], end - start)
                ax = genes.call_args.args[0]
                self.assertIn('GENE1', [text.get_text() for text in ax.texts])
                title_chrom = ax.get_title().split()[1][:-1]
                self.assertEqual(title_chrom, graph_chrom)
                if not with_cnr:
                    self.assertNotIn('--cnr None', (self.root / 'plot_visualization.log').read_text())

    def test_custom_centromere_is_used_with_either_prefix_style(self):
        for graph_chrom, bed_chrom in [('7', 'chr7'), ('chr7', '7')]:
            with self.subTest(graph=graph_chrom, bed=bed_chrom):
                graph, cycles, _, _, _, _ = self.fixtures(graph_chrom)
                bed = self.root / 'centromere.bed'
                bed.write_text(f'{bed_chrom}\t500000\t500100\n')
                with patch.object(viz.plt, 'title', wraps=viz.plt.title) as title:
                    viz.visualize_BFB(cycles, graph, None, str(self.root / 'plot'), centromere=bed)
                self.assertEqual(title.call_args.args[0], f'plot {graph_chrom}q')

    def test_empty_cnr_match_reports_the_chromosome(self):
        graph, cycles, cnr, _, _, _ = self.fixtures('7', '8')
        with self.assertRaisesRegex(ValueError, "No copy-number bins for chromosome '7'"):
            viz.visualize_BFB(cycles, graph, cnr, str(self.root / 'plot'))


if __name__ == '__main__':
    unittest.main()
