import unittest
import os
import logging
import filecmp
from albatradis.ScatterPlot import ScatterPlot

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','scatterplot')


class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

class TestScatterPlot(unittest.TestCase):

	def test_valid(self):
		s = ScatterPlot([os.path.join(data_dir, 'condition1.plot'), os.path.join(data_dir, 'condition2.plot')], [os.path.join(data_dir, 'control1.plot'), os.path.join(data_dir, 'control2.plot')], 10,  'scattertest', True, verbose = True)
		self.assertTrue(s.create_scatter_plot())
		self.assertTrue(s.create_linear_plot())
		self.assertTrue(s.create_abs_scatter_plot())

