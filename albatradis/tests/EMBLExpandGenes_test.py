import unittest
import os
import filecmp
import logging
from albatradis.EMBLExpandGenes import EMBLExpandGenes

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','emblexpandgenes')

class TestEMBLExpandGenes(unittest.TestCase):

	def test_embl_expand_genes(self):
		e = EMBLExpandGenes(os.path.join(data_dir, 'one_gene') , 5)
		
		output_file = os.path.join(data_dir, 'actual_one_gene')
		e.construct_file(output_file)
		self.assertTrue(os.path.exists(output_file))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_one_gene'), output_file))
		
		os.remove(output_file)	
		
		