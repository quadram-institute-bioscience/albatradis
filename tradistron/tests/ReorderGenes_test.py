import unittest
import os
import logging
import filecmp
from tradistron.ReorderGenes import ReorderGenes
from tradistron.GeneToFiles import GeneToFiles

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

class TestNormalisePlots(unittest.TestCase):

	def test_valid(self):
		r = ReorderGenes(['geneA','geneB','geneC'], ['file1','file2','file3','file4'], {'file1': [0,0,9],'file2': [9,-9,9], 'file3': [0,5,9], 'file4': [1, 0, 1]})
		self.assertEqual('geneC', r.get_highest_freq(r.genes_to_files))
		
	def test_zeros(self):
		r = ReorderGenes(['geneA','geneB','geneC'], ['file1','file2','file3','file4'], {'file1': [0,0,0],'file2': [0,0,0], 'file3': [0,0,0], 'file4': [0, 0, 0]})
		self.assertEqual('geneA', r.get_highest_freq(r.genes_to_files))

	def test_common_counting(self):
		g1 = GeneToFiles('abc', gene_to_files= [0,11,13,0])
		g2 = GeneToFiles('efg', gene_to_files= [0,11,-12,0])
		
		r = ReorderGenes([], [], {})
		self.assertEqual(2,r.files_in_common(g1, g2))
		
	def test_find_closest_gene(self):
		r = ReorderGenes([], [], {})
		gf = {'abc': GeneToFiles('abc', gene_to_files= [0,2,2,0]), 'efg': GeneToFiles('efg', gene_to_files= [2,2,-2,0]), 'geneA': GeneToFiles('geneA', gene_to_files= [0,2,2,0]),  'outlier': GeneToFiles('outlier', gene_to_files= [0,0,0,1]) }
		
		self.assertEqual('abc', r.find_closest_gene( 'geneA', gf))
		self.assertEqual('abc', r.find_closest_gene( 'efg', gf))
		self.assertEqual('efg', r.find_closest_gene( 'abc', gf))
		self.assertEqual('abc', r.find_closest_gene( 'outlier', gf))
		
	def test_reorder_genes(self):
		r = ReorderGenes(['geneA','geneB','geneC'], ['file1','file2','file3','file4'], {'file1': [0,0,9],'file2': [9,-9,9], 'file3': [0,5,9], 'file4': [1, 0, 1]})
		self.assertEqual(['geneC', 'geneA', 'geneB'], r.reorder_genes())
		