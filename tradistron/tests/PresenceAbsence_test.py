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
		
		if not os.path.exists('testoutput' ):
			os.makedirs('testoutput' )	
		
		genereports = [os.path.join(data_dir, 'ctrl.csv'), os.path.join(data_dir, 'gent1.csv'), os.path.join(data_dir, 'gent5.csv'), os.path.join(data_dir, 'gent25.csv') ]
		emblfile = os.path.join(data_dir, 'reference.embl')
		
		all_outputfile = os.path.join('testoutput', 'all_logfc.csv')
		filtered_outputfile =os.path.join('testoutput', 'filtered_logfc.csv')
		p = PresenceAbsence(genereports, emblfile,  False, 'testoutput')
		
		self.assertEqual(962, len(p.gene_names))
		self.assertEqual(1079, len(p.features))
		self.assertTrue(p.create_output_files())
		
		
		self.assertTrue(all_outputfile)
		self.assertTrue(filtered_outputfile)
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_logfc.csv'), all_outputfile))
		
		shutil.rmtree('testoutput')
		