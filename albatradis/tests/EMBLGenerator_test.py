import unittest
import os
import logging
from albatradis.EMBLGenerator import EMBLGenerator
from albatradis.WindowGenerator import WindowGenerator

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','emblgenerator')

class TestEMBLGenerator(unittest.TestCase):

	def test_embl_generator_nonoverlap(self):
		w = WindowGenerator(20, 5,5)
		e = EMBLGenerator(w.create_windows(), 20)
		
		e.construct_file(os.path.join(data_dir, 'nonoverlap'))
		self.assertTrue(os.path.exists(os.path.join(data_dir, 'nonoverlap')))
		
