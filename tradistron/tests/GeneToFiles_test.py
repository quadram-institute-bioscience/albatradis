import unittest
import os
import logging
import filecmp
from tradistron.GeneToFiles import GeneToFiles

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

class TestNormalisePlots(unittest.TestCase):

	def test_empty(self):
		g = GeneToFiles('abc', gene_to_files= [0,0,0,0])
		self.assertEqual(0, g.number_of_files)
		
	def test_one(self):
		g = GeneToFiles('abc', gene_to_files= [0,11,0,0])
		self.assertEqual(1, g.number_of_files)
		
	def test_all(self):
		g = GeneToFiles('abc', gene_to_files= [2,11,10,-15])
		self.assertEqual(4, g.number_of_files)
		
