import unittest
import os
import logging
import filecmp
from albatradis.HeatMap import HeatMap

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','heatmap')

class TestNormalisePlots(unittest.TestCase):

	def test_same_path(self):
		output_file = os.path.join(data_dir,'heatmap.png')
		reports_to_logfc = {'/path/Antibiotic_1.csv': [1,5,0], '/path/Antibiotic_2': [0,10,-10], '/path/someotherAntibiotic_0.5': [0,1,0]}
		h = HeatMap(reports_to_logfc, ['geneA','geneB', 'geneC'], output_file)


		self.assertEqual(['Antibiotic 1', 'Antibiotic 2', 'someotherAntibiotic 0.5'], h.clean_filenames())
		self.assertTrue(h.create_heat_map())
		self.assertTrue(os.path.exists(output_file))
		os.remove(output_file)
		
	def test_no_common_path(self):
		output_file = os.path.join(data_dir,'heatmap.png')
		reports_to_logfc = {'/path123/Antibiotic_1.csv': [1,5,0], '/path456/Antibiotic_2': [0,10,-10], '/path/someotherAntibiotic_0.5': [0,1,0]}
		h = HeatMap(reports_to_logfc, ['geneA','geneB', 'geneC'], output_file)

		self.assertEqual(['/path123/Antibiotic 1', '/path456/Antibiotic 2', '/path/someotherAntibiotic 0.5'], h.clean_filenames())
		self.assertTrue(h.create_heat_map())
		self.assertTrue(os.path.exists(output_file))
		os.remove(output_file)
		
