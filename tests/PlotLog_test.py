import os
import unittest

from albatradis.PlotLog import PlotLog

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'plotlog')


class TestPlotLog(unittest.TestCase):

    def setUp(self) -> None:
        self.p = PlotLog('x',
                         genome_length=20,
                         minimum_logfc=2,
                         maximum_pvalue=0.05,
                         maximum_qvalue=0.05,
                         minimum_logcpm=1,
                         window_size=4,
                         span_gaps=1,
                         report_decreased_insertions=True,
                         embl_file=None)

    def test_span_gaps_merge(self):
        l = self.p.span_block_gaps([0, 0, 0, 9, 9, 9, 9, 0, 0, 0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 0],
                                   [0.0] * 20,
                                   [0.0] * 20)
        self.assertEqual([0, 0, 0, 9, 9, 9, 9, 2, 2, 2, 2, 8, 8, 8, 8, 0, 0, 0, 0, 0], l[0])

    def test_span_gaps_nomerge(self):
        l = self.p.span_block_gaps([0, 0, 0, 9, 9, 9, 0, 0, 0, 0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 0],
                                   [0.0] * 20,
                                   [0.0] * 20)
        self.assertEqual([0, 0, 0, 9, 9, 9, 0, 0, 0, 0, 0, 8, 8, 8, 8, 0, 0, 0, 0, 0], l[0])

    def test_span_gaps_block_at_end(self):
        l = self.p.span_block_gaps([0, 0, 0, 9, 9, 9, 0, 0, 0, 0, 0, 8, 8, 8, 8, 0, 0, 7, 7, 7],
                                   [0.0] * 20,
                                   [0.0] * 20)
        self.assertEqual([0, 0, 0, 9, 9, 9, 0, 0, 0, 0, 0, 8, 8, 8, 8, 2, 2, 7, 7, 7], l[0])

    def test_span_gaps_merge_neg(self):
        l = self.p.span_block_gaps([0, 0, 0, -9, -9, -9, -9, 0, 0, 0, 0, -8, -8, -8, -8, 0, 0, 0, 0, 0],
                                   [0.0] * 20,
                                   [0.0] * 20)
        self.assertEqual([0, 0, 0, -9, -9, -9, -9, -2, -2, -2, -2, -8, -8, -8, -8, 0, 0, 0, 0, 0], l[0])

    def test_span_gaps_nomerge_neg(self):
        l = self.p.span_block_gaps([0, 0, 0, -9, -9, -9, 0, 0, 0, 0, 0, -8, -8, -8, -8, 0, 0, 0, 0, 0],
                                   [0.0] * 20,
                                   [0.0] * 20)
        self.assertEqual([0, 0, 0, -9, -9, -9, 0, 0, 0, 0, 0, -8, -8, -8, -8, 0, 0, 0, 0, 0], l[0])

    def test_span_gaps_block_at_end_neg(self):
        l = self.p.span_block_gaps([0, 0, 0, -9, -9, -9, 0, 0, 0, 0, 0, -8, -8, -8, -8, 0, 0, -7, -7, -7],
                                   [0.0] * 20,
                                   [0.0] * 20)
        self.assertEqual([0, 0, 0, -9, -9, -9, 0, 0, 0, 0, 0, -8, -8, -8, -8, -2, -2, -7, -7, -7], l[0])
