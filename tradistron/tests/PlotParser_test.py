import unittest
import os
import logging
from tradistron.PlotParser import PlotParser

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','plotparser')

class TestPlotParser(unittest.TestCase):

	def test_valid_unziped(self):
		d = PlotParser(os.path.join(data_dir,'valid'), 0)
		
		self.assertEqual(d.forward, [0,0,0,0,0,1,1,3])
		self.assertEqual(d.reverse, [0,1,4,5,0,0,0,0])
		self.assertEqual(d.combined, [0,1,4,5,0,1,1,3])
		
	def test_valid_ziped(self):
		d = PlotParser(os.path.join(data_dir,'valid.gz'),0)

		self.assertEqual(d.forward, [0,0,0,0,0,1,1,3])
		self.assertEqual(d.reverse, [0,1,4,5,0,0,0,0])
		self.assertEqual(d.combined, [0,1,4,5,0,1,1,3])

	def test_valid_threshold(self):
		d = PlotParser(os.path.join(data_dir,'valid'), 3)

		self.assertEqual(d.forward, [0,0,0,0,0,0,0,3])
		self.assertEqual(d.reverse, [0,0,4,5,0,0,0,0])
		self.assertEqual(d.combined, [0,0,4,5,0,0,0,3])