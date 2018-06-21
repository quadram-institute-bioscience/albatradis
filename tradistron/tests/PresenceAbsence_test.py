import unittest
import os
import filecmp
import shutil
from tradistron.PresenceAbsence import PresenceAbsence

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

data_dir = os.path.join('tradistron','tests', 'data','presenceabsence')

class TestPresenceAbsence(unittest.TestCase):

	def test_valid_file(self):
		genereports = [os.path.join(data_dir, 'ctrl.csv'), os.path.join(data_dir, 'gent1.csv'), os.path.join(data_dir, 'gent5.csv'), os.path.join(data_dir, 'gent25.csv') ]
		emblfile = os.path.join(data_dir, 'reference.embl')
		
		outputfile = os.path.join(data_dir, 'logfc.csv')
		p = PresenceAbsence(genereports, emblfile, False, False)
		
		self.assertEqual(962, len(p.gene_names))
		self.assertEqual(1079, len(p.features))
		self.assertTrue(p.create_gene_logfc_spreadsheet(outputfile))
		self.assertTrue(os.path.exists(outputfile))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_logfc.csv'), outputfile))
		
		os.remove(outputfile)
		