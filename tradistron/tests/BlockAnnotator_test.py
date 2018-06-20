import unittest
import os
import logging
import filecmp
from tradistron.BlockAnnotator import BlockAnnotator
from tradistron.Block import Block

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','blockannotator')

class TestBlockAnnotator(unittest.TestCase):

	def test_fully_contained(self):
		embl_file = os.path.join(data_dir,'annotation.embl')
		blocks = [Block(1,110, 110, 10, 'x')]
		a = BlockAnnotator(embl_file, blocks)
		a.annotate_blocks()
		self.assertEqual(1, len(blocks[0].genes))
		self.assertEqual('total_inactivation', blocks[0].genes[0].category)
		
	def test_fully_contained_two(self):
		embl_file = os.path.join(data_dir,'annotation.embl')
		blocks = [Block(1,200, 110, 10, 'x')]
		a = BlockAnnotator(embl_file, blocks)
		a.annotate_blocks()
		self.assertEqual(2, len(blocks[0].genes))
		self.assertEqual('total_inactivation', blocks[0].genes[0].category)

	def test_total_file_in_block(self):
		embl_file = os.path.join(data_dir,'annotation.embl')
		blocks = [Block(1,600, 600, 10, 'x')]
		a = BlockAnnotator(embl_file, blocks)
		a.annotate_blocks()
		self.assertEqual(4, len(blocks[0].genes))
		self.assertEqual('total_inactivation', blocks[0].genes[0].category)

	def test_block_near_end(self):
		embl_file = os.path.join(data_dir,'annotation.embl')
		blocks = [Block(90,110, 20, 10, 'x')]
		blocks[0].direction = 'reverse'
		a = BlockAnnotator(embl_file, blocks)
		a.annotate_blocks()
		self.assertEqual(1, len(blocks[0].genes))
		self.assertEqual('increased_mutants_at_end_of_gene', blocks[0].genes[0].category)
		
	def test_block_near_end_reverse(self):
		embl_file = os.path.join(data_dir,'reverse.embl')
		blocks = [Block(1,20, 20, 10, 'x')]
		blocks[0].direction = 'forward'
		a = BlockAnnotator(embl_file, blocks)
		a.annotate_blocks()
		self.assertEqual(1, len(blocks[0].genes))
		self.assertEqual('increased_mutants_at_end_of_gene', blocks[0].genes[0].category)
		
	def test_block_near_end_decrease(self):
		embl_file = os.path.join(data_dir,'annotation.embl')
		blocks = [Block(90,110, 20, -10, 'x')]
		blocks[0].direction = 'reverse'
		a = BlockAnnotator(embl_file, blocks)
		a.annotate_blocks()
		self.assertEqual(1, len(blocks[0].genes))
		self.assertEqual('decreased_mutants_at_end_of_gene', blocks[0].genes[0].category)
		self.assertEqual(['1_100'], blocks[0].genes[0].feature.qualifiers["gene"])
	