import unittest
import os
import logging
from tradistron.GeneReport import GeneReport

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','genereport')

class TestGeneReport(unittest.TestCase):

	def test_valid_file(self):
		g = GeneReport(os.path.join(data_dir, 'ctrl.csv'))
		self.assertTrue(g)
		gene_names  = ['perR','ykfC','ykgH','yaiL','yaiP','unknown']
		self.assertEqual(['10', '11', '10', '-10', '11', '10'], g.genes_to_logfc(gene_names))
		
		# reorder the gene names
		gene_names  = ['yaiL','ykfC','ykgH','perR','unknown', 'yaiP']
		self.assertEqual(['-10', '11', '10', '10', '10', '11'], g.genes_to_logfc(gene_names))
		