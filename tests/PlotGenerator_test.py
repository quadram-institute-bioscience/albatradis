import filecmp
import os
import unittest

from albatradis.PlotGenerator import PlotGenerator
from albatradis.PlotParser import PlotParser

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'plotgenerator')

class TestPlotGenerator(unittest.TestCase):

	def test_forward_only(self):
		filename = os.path.join(data_dir, 'test_plotgen')
		p = PlotGenerator([0,0,0,0,0,1,1,3], [], filename)
		self.assertTrue(p.construct_file())
		
		self.assertTrue(os.path.exists(filename))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_plotgen'), filename))
		os.remove(filename)
		
	def test_generate_and_read(self):
		filename = os.path.join(data_dir, 'plot.test')
		p = PlotGenerator([0,0,0,0,0,1,1,3], [9,9,0,9,9,1,1,3], filename)
		p.construct_file()
		
		d = PlotParser(filename,0)
		self.assertTrue(self.check_arrays_equal(d.forward, [0,0,0,0,0,1,1,3]))
		self.assertTrue(self.check_arrays_equal(d.reverse, [9,9,0,9,9,1,1,3]))
		self.assertTrue(self.check_arrays_equal(d.combined, [9,9,0,9,9,2,2,6]))
		
		os.remove(filename)
		
	def check_arrays_equal(self, array1, array2):
		for i, val in enumerate(array1):
			if array1[i] != array2[i]:
				print(str(array1[i])+" not equal to "+str(array2[i]))
				return False
		return True