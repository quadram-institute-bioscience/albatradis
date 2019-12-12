import unittest
import os
import filecmp
import shutil
import subprocess


from albatradis.PresenceAbsence import PresenceAbsence

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join('data','presenceabsence')

base_dir = os.path.abspath(os.path.join(test_modules_dir, '..', '..'))

example_dir = os.path.join('/albatradis', 'data', 'presence_absence_data')


class TestPresenceAbsence(unittest.TestCase):

	def test_valid_file(self):
		
		if not os.path.exists('testoutput' ):
			os.makedirs('testoutput' )	
		
		genereports = [os.path.join(data_dir, 'ctrl.csv'), os.path.join(data_dir, 'gent1.csv'), os.path.join(data_dir, 'gent5.csv'), os.path.join(data_dir, 'gent25.csv') ]
		emblfile = os.path.join(data_dir, 'reference.embl')

		all_outputfile = os.path.join('testoutput', 'all_logfc.csv')
		filtered_outputfile = os.path.join('testoutput', 'filtered_logfc.csv')
		dendrogram = os.path.join('testoutput', 'distance_matrix_dendrogram.tre')
		nj_tree = os.path.join('testoutput', 'nj_newick_tree.tre')
		p = PresenceAbsence(genereports, emblfile,  False, 'testoutput')
		exp_lfc=os.path.join(data_dir, 'expected_logfc.csv')

		self.assertEqual(962, len(p.gene_names))
		self.assertEqual(1079, len(p.features))
		self.assertTrue(p.create_output_files())
		self.assertTrue(os.path.exists(exp_lfc))
		
		self.assertTrue(all_outputfile)
		self.assertTrue(filtered_outputfile)
		self.assertTrue(dendrogram)
		self.assertTrue(nj_tree)
		self.assertTrue(filecmp.cmp(exp_lfc, all_outputfile))
		
		#sshutil.rmtree('testoutput')
		
	def test_example_toy(self):
		files = " ".join([os.path.join(example_dir, 'gene_report_1mgL.csv'), os.path.join(example_dir, 'gene_report_003mgL.csv'), os.path.join(example_dir, 'gene_report_05mgL.csv'), os.path.join(example_dir, 'gene_report_006mgL.csv'), os.path.join(example_dir, 'gene_report_0008mgL.csv'), os.path.join(example_dir, 'gene_report_0015mgL.csv'), os.path.join(example_dir, 'gene_report_025mgL.csv'), os.path.join(example_dir, 'gene_report_0125mgL.csv')])
		emblfile = os.path.join(example_dir, 'reference_BW25113.embl')
		out_dir = (os.path.join(base_dir, 'albatradis_output:/work'))

		cmd = " ".join(
			['docker run --rm  -v', out_dir, 'quadraminstitute/albatradis:latest albatradis-presence_absence ', emblfile, files])

		subprocess.call(cmd, shell=True)

		self.assertTrue(os.path.exists(os.path.join(base_dir, 'albatradis_output/output')))
		self.assertTrue(os.path.exists(os.path.join(base_dir, 'albatradis_output/output/logfc.dot')))
		shutil.rmtree(os.path.join(base_dir,'albatradis_output'))