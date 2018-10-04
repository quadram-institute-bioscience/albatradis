import unittest
import os
import logging
import filecmp
from albatradis.PlotLog import PlotLog

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','plotlog')

class TestPlotLog(unittest.TestCase):
	
	def test_span_gaps_merge(self):
		p = PlotLog('x', 20, 2, 0.05, 1, 4, 1, True, None)
		l = p.span_block_gaps([0,0,0,9,9,9,9,0,0,0,0,8,8,8,8,0,0,0,0,0])
		self.assertEqual([0,0,0,9,9,9,9,2,2,2,2,8,8,8,8,0,0,0,0,0],l )

	def test_span_gaps_nomerge(self):
		p = PlotLog('x', 20, 2, 0.05, 1, 4, 1, True, None)
		l = p.span_block_gaps([0,0,0,9,9,9,0,0,0,0,0,8,8,8,8,0,0,0,0,0])
		self.assertEqual([0,0,0,9,9,9,0,0,0,0,0,8,8,8,8,0,0,0,0,0],l )

	def test_span_gaps_block_at_end(self):
		p = PlotLog('x', 20, 2, 0.05, 1, 4, 1, True, None)
		l = p.span_block_gaps([0,0,0,9,9,9,0,0,0,0,0,8,8,8,8,0,0,7,7,7])
		self.assertEqual([0,0,0,9,9,9,0,0,0,0,0,8,8,8,8,2,2,7,7,7],l )
		
	def test_span_gaps_merge_neg(self):
		p = PlotLog('x', 20, 2, 0.05, 1, 4, 1, True, None)
		l = p.span_block_gaps([0,0,0,-9,-9,-9,-9,0,0,0,0,-8,-8,-8,-8,0,0,0,0,0])
		self.assertEqual([0,0,0,-9,-9,-9,-9,-2,-2,-2,-2,-8,-8,-8,-8,0,0,0,0,0],l )

	def test_span_gaps_nomerge_neg(self):
		p = PlotLog('x', 20, 2, 0.05, 1, 4, 1, True, None)
		l = p.span_block_gaps([0,0,0,-9,-9,-9,0,0,0,0,0,-8,-8,-8,-8,0,0,0,0,0])
		self.assertEqual([0,0,0,-9,-9,-9,0,0,0,0,0,-8,-8,-8,-8,0,0,0,0,0],l )

	def test_span_gaps_block_at_end_neg(self):
		p = PlotLog('x', 20, 2, 0.05, 1, 4, 1, True, None)
		l = p.span_block_gaps([0,0,0,-9,-9,-9,0,0,0,0,0,-8,-8,-8,-8,0,0,-7,-7,-7])
		self.assertEqual([0,0,0,-9,-9,-9,0,0,0,0,0,-8,-8,-8,-8,-2,-2,-7,-7,-7],l )
	