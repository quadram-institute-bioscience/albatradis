import unittest
import os
import logging
import filecmp
from albatradis.ArtemisProject import ArtemisProject

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

data_dir = os.path.join('albatradis','tests', 'data','experimentcollection')

class TestArtemisProject(unittest.TestCase):

	def test_vary_mic(self):
		outputfile = "project_file"
		ap = ArtemisProject(outputfile, False, os.path.join(data_dir, 'spreadsheet'), os.path.join(data_dir, 'reference'), [os.path.join(data_dir, 'control')])
		
		self.assertTrue(ap.create_project_file())
		self.assertTrue(os.path.exists(outputfile))
		
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_project_file'), outputfile))
		os.remove(outputfile)
