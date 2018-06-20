import unittest
import os
import logging
import filecmp
from tradistron.GeneAnnotator import GeneAnnotator
from tradistron.Block import Block

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','geneannotator')

class TestGeneAnnotator(unittest.TestCase):

	def test_fully_contained(self):
		embl_file = os.path.join(data_dir,'annotation.embl')
		blocks = [Block(1,110, 110, 10, 'x')]
		a = GeneAnnotator(embl_file, blocks)
		genes = a.annotate_genes()
		self.assertEqual(1, len(genes))
		self.assertEqual('total_inactivation', genes[0].category())
		
	def test_fully_contained_three(self):
		embl_file = os.path.join(data_dir,'annotation.embl')
		blocks = [Block(1,200, 110, 10, 'x')]
		a = GeneAnnotator(embl_file, blocks)
		genes = a.annotate_genes()
		self.assertEqual(3, len(genes))
		self.assertEqual('total_inactivation', genes[0].category())

	def test_total_file_in_block(self):
		embl_file = os.path.join(data_dir,'annotation.embl')
		blocks = [Block(1,600, 600, 10, 'x')]
		a = GeneAnnotator(embl_file, blocks)
		genes = a.annotate_genes()
		self.assertEqual(4, len(genes))
		self.assertEqual('total_inactivation', genes[0].category())

	def test_block_near_end(self):
		embl_file = os.path.join(data_dir,'annotation.embl')
		blocks = [Block(90,110, 20, 10, 'x')]
		blocks[0].direction = 'reverse'
		a = GeneAnnotator(embl_file, blocks)
		genes = a.annotate_genes()
		self.assertEqual(1, len(genes))
		self.assertEqual('increased_mutants_at_end_of_gene', genes[0].category())
		
	def test_block_near_end_reverse(self):
		embl_file = os.path.join(data_dir,'reverse.embl')
		blocks = [Block(1,20, 20, 10, 'x')]
		blocks[0].direction = 'forward'
		a = GeneAnnotator(embl_file, blocks)
		genes = a.annotate_genes()
		self.assertEqual(1, len(genes))
		self.assertEqual('increased_mutants_at_end_of_gene', genes[0].category())
		
	def test_block_near_end_decrease(self):
		embl_file = os.path.join(data_dir,'annotation.embl')
		blocks = [Block(90,110, 20, -10, 'x')]
		blocks[0].direction = 'reverse'
		a = GeneAnnotator(embl_file, blocks)
		genes = a.annotate_genes()
		self.assertEqual(1, len(genes))
		self.assertEqual('decreased_mutants_at_end_of_gene', genes[0].category())
		self.assertEqual('1_100', genes[0].gene_name())
	