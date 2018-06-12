import unittest
import os
import logging
from tradistron.TradisTron import TradisTron
import shutil


test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','tradistron')

class TestOptions:
	def __init__(self, plotfiles, minimum_threshold, window_size,window_interval, verbose, prefix, minimum_logcpm, minimum_logfc, pvalue, iterations ):
		self.plotfiles = plotfiles
		self.minimum_threshold = minimum_threshold
		self.window_size = window_size
		self.window_interval = window_interval
		self.verbose = verbose	
		self.prefix = prefix
		self.minimum_logcpm = minimum_logcpm
		self.minimum_logfc = minimum_logfc
		self.pvalue = pvalue
		self.iterations = iterations

class TestTradisTron(unittest.TestCase):
	
	def test_small_real(self):
		case = os.path.join(data_dir, 'small_case.insert_site_plot.gz')
		control = os.path.join(data_dir, 'small_control.insert_site_plot.gz')
		t = TradisTron(TestOptions([case, control], 3, 100, 100, False, 'testoutput', 1, 1, 1, 1))
		self.assertTrue(t.run())
		self.assertTrue(os.path.exists('testoutput_1'))
		shutil.rmtree("testoutput_1")
		
		
	def test_small_2iterations(self):
		case = os.path.join(data_dir, 'small_case.insert_site_plot.gz')
		control = os.path.join(data_dir, 'small_control.insert_site_plot.gz')
		t = TradisTron(TestOptions([case, control], 3, 100, 100, False, 'testoutput', 1, 1, 1, 2))
		self.assertTrue(t.run())
		self.assertTrue(os.path.exists('testoutput_1'))
		shutil.rmtree("testoutput_1")
		self.assertTrue(os.path.exists('testoutput_2'))
		shutil.rmtree("testoutput_2")
		
	#def test_big_real(self):
	#	case = os.path.join(data_dir, 'big_case.insert_site_plot.gz')
	#	control = os.path.join(data_dir, 'big_control.insert_site_plot.gz')
	#	t = TradisTron(TestOptions([case, control], 3, 100, 100, False, 'testoutput', 1, 1,1,1 ))
	#	self.assertTrue(t.run())
	#	self.assertTrue(os.path.exists('testoutput_1'))
	#	shutil.rmtree("testoutput_1")