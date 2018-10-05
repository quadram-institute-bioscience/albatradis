import unittest
import os
import logging
import filecmp
from albatradis.PrepareEMBLFile import PrepareEMBLFile

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','prepareinputfiles')
data_expand_genes_dir = os.path.join(test_modules_dir, 'data','emblexpandgenes')

class TestPrepareEMBLFile(unittest.TestCase):

	def test_small_valid_file(self):
		p = PrepareEMBLFile(os.path.join(data_dir,'valid'), 1, 4, 2, False, 100, None)
		
		self.assertTrue(p.create_file())
		self.assertTrue(os.path.exists(p.embl_filename))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_embl_filename'), p.embl_filename))

		os.remove(p.embl_filename)	
		
	def test_reuse_annotation(self):
		p = PrepareEMBLFile(os.path.join(data_dir,'valid'), 1, 4, 2, True, 15, os.path.join(data_expand_genes_dir, 'one_gene'))
		self.assertTrue(p.create_file())
		
		self.assertTrue(filecmp.cmp(os.path.join(data_expand_genes_dir, 'expected_one_gene'), p.embl_filename))
		os.remove(p.embl_filename)
