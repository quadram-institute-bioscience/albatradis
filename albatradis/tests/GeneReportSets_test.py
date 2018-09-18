import unittest
import os
import logging
import filecmp
import shutil
from tradistron.GeneReportSets import GeneReportSets

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','genereportsets')

class TestGeneReportSets(unittest.TestCase):

	def test_same_file_twice(self):
		g = GeneReportSets([os.path.join(data_dir, 'sample1.csv'), os.path.join(data_dir, 'sample1.csv')], os.path.join(data_dir, 'same_file'))
		self.assertTrue(g)
		self.assertTrue(g.write_union_file())
		self.assertTrue(os.path.exists(os.path.join(data_dir, 'same_file','union_gene_report.csv')))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'same_file' ,'union_gene_report.csv'), os.path.join(data_dir, 'expected_same_file_' +'union_gene_report.csv')))
		shutil.rmtree(os.path.join(data_dir, 'same_file'))
	
	def test_two_files_union(self):
		g = GeneReportSets([os.path.join(data_dir, 'sample1.csv'), os.path.join(data_dir, 'sample2.csv')], os.path.join(data_dir, 'two_files'))
		self.assertTrue(g)
		self.assertTrue(g.write_union_file())
		self.assertTrue(os.path.exists(os.path.join(data_dir, 'two_files','union_gene_report.csv')))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'two_files' ,'union_gene_report.csv'), os.path.join(data_dir, 'expected_two_files_' +'union_gene_report.csv')))
		shutil.rmtree(os.path.join(data_dir, 'two_files'))

	def test_two_files_intersection(self):
		g = GeneReportSets([os.path.join(data_dir, 'sample1.csv'), os.path.join(data_dir, 'sample2.csv')], os.path.join(data_dir, 'two_filesi'))
		self.assertTrue(g)
		self.assertTrue(g.write_intersection_file())
		print(os.path.join(data_dir, 'two_filesi'))
		self.assertTrue(os.path.exists(os.path.join(data_dir, 'two_filesi','intersection_gene_report.csv')))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'two_filesi' ,'intersection_gene_report.csv'), os.path.join(data_dir, 'expected_two_files_' +'intersection_gene_report.csv')))
		
		shutil.rmtree(os.path.join(data_dir, 'two_filesi'))