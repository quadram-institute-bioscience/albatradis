import unittest
import os
import logging
import filecmp
from tradistron.PlotMasking import PlotMasking

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','plotmasking')

class TestPlotMasking(unittest.TestCase):
	
	def test_plot_masking_small(self):
		case = os.path.join(data_dir, 'sample1.plot')
		control = os.path.join(data_dir, 'sample2.plot')
		mask = os.path.join(data_dir, 'mask.plot')
		pm = PlotMasking([case, control], mask)
		self.assertEqual(2, len(pm.output_plot_files))
		
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_sample1.plot'), pm.output_plot_files[case]))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_sample2.plot'), pm.output_plot_files[control]))
