import unittest
import os
import logging
import filecmp
from albatradis.PlotGenerator import PlotGenerator

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','plotgenerator')

class TestPlotGenerator(unittest.TestCase):

	def test_forward_only(self):
		filename = os.path.join(data_dir, 'test_plotgen')
		p = PlotGenerator([0,0,0,0,0,1,1,3], [], filename)
		self.assertTrue(p.construct_file())
		
		self.assertTrue(os.path.exists(filename))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_plotgen'), filename))
		os.remove(filename)
