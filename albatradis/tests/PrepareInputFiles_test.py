import unittest
import os
import logging
import filecmp
from albatradis.PrepareInputFiles import PrepareInputFiles

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','prepareinputfiles')
data_expand_genes_dir = os.path.join(test_modules_dir, 'data','emblexpandgenes')

class TestPrepareInputFiles(unittest.TestCase):

	def test_small_valid_file(self):
		p = PrepareInputFiles(os.path.join(data_dir,'valid'), 1)
		
		self.assertTrue(p.create_all_files())
		self.assertTrue(os.path.exists(p.forward_plot_filename))
		self.assertTrue(os.path.exists(p.reverse_plot_filename))
		self.assertTrue(os.path.exists(p.combined_plot_filename))
		
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_forward_plot_filename'), p.forward_plot_filename))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_reverse_plot_filename'), p.reverse_plot_filename))
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_combined_plot_filename'), p.combined_plot_filename))
		
		os.remove(p.forward_plot_filename)
		os.remove(p.reverse_plot_filename)
		os.remove(p.combined_plot_filename)
