import unittest
import os
import logging
from albatradis.WindowGenerator import WindowGenerator

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data', 'windowgenerator')

class TestWindowGenerator(unittest.TestCase):

	def test_window_generator_nonoverlap(self):
		w = WindowGenerator(20, 5,5)
		self.assertEqual(4, len(w.create_windows()))
		
	def test_window_generator_overlap(self):
		w = WindowGenerator(20, 5,1)
		self.assertEqual(16, len(w.create_windows()))
	
