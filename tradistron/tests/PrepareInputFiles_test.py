import unittest
import os
import logging
from tradistron.PrepareInputFiles import PrepareInputFiles

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','prepareinputfiles')

class TestPrepareInputFiles(unittest.TestCase):

	def test_small_valid_file(self):
		p = PrepareInputFiles(os.path.join(data_dir,'valid'), 1, 4, 2)
		
		self.assertTrue(p.create_all_files())
