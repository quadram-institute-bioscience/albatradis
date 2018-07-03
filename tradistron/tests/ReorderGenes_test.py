import unittest
import os
import logging
import filecmp
from tradistron.ReorderGenes import ReorderGenes

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

class TestNormalisePlots(unittest.TestCase):

	def test_valid(self):
		r = ReorderGenes(['geneA','geneB','geneC'], ['file1','file2','file3','file4'], {'file1': [0,0,9],'file2': [9,-9,9], 'file3': [0,5,9], 'file4': [1, 0, 1]})
		self.assertEqual(4, r.get_highest_freq(r.genes_to_files))
		
	def test_zeros(self):
		r = ReorderGenes(['geneA','geneB','geneC'], ['file1','file2','file3','file4'], {'file1': [0,0,0],'file2': [0,0,0], 'file3': [0,0,0], 'file4': [0, 0, 0]})
		self.assertEqual(0, r.get_highest_freq(r.genes_to_files))
