import unittest
import os
import logging
from tradistron.TradisTron import TradisTron

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','tradistron')

class TestOptions:
	def __init__(self, plotfiles, minimum_threshold, window_size,window_interval, verbose ):
		self.plotfiles = plotfiles
		self.minimum_threshold = minimum_threshold
		self.window_size = window_size
		self.window_interval = window_interval
		self.verbose = verbose	

class TestTradisTron(unittest.TestCase):
	
	def test_small_real(self):
		case = os.path.join(data_dir, 'small_case.insert_site_plot.gz')
		control = os.path.join(data_dir, 'small_control.insert_site_plot.gz')
		t = TradisTron(TestOptions([case, control], 3, 100, 100, False))
		self.assertTrue(t.run())
		