import unittest
import os
import logging
import filecmp
from tradistron.NormalisePlots import NormalisePlots

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','normaliseplots')

class TestNormalisePlots(unittest.TestCase):

	def test_big_differences(self):
		p = NormalisePlots([os.path.join(data_dir,'sample1'), os.path.join(data_dir,'sample2')])
		output_files = p.create_normalised_files()
		self.assertEqual(2, len(output_files))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'sample2'), output_files[1]))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_sample1'), output_files[0]))
		