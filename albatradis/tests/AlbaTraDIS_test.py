import filecmp
import os
import shutil
import unittest
import subprocess



from albatradis.AlbaTraDIS import AlbaTraDIS

test_modules_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(test_modules_dir, 'data','albatradis')
base_dir = os.path.abspath(os.path.join(test_modules_dir, '..', '..'))
example_dir = os.path.join('/albatradis', 'data', 'albatradis_data')

class TestOptions:
	def __init__(self, plotfiles, minimum_threshold, window_size,window_interval, verbose, prefix, minimum_logcpm, minimum_logfc, pvalue, qvalue, iterations, dont_normalise_plots,minimum_block,span_gaps, emblfile, minimum_proportion_insertions, strict_signal, use_annotation, prime_feature_size):
		self.plotfiles = plotfiles
		self.minimum_threshold = minimum_threshold
		self.window_size = window_size
		self.window_interval = window_interval
		self.verbose = verbose	
		self.prefix = prefix
		self.minimum_logcpm = minimum_logcpm
		self.minimum_logfc = minimum_logfc
		self.pvalue = pvalue
		self.qvalue = qvalue
		self.iterations = iterations
		self.dont_normalise_plots = dont_normalise_plots
		self.minimum_block = minimum_block
		self.span_gaps = span_gaps
		self.emblfile = emblfile
		self.minimum_proportion_insertions = minimum_proportion_insertions
		self.strict_signal = strict_signal
		self.use_annotation    = use_annotation
		self.prime_feature_size = prime_feature_size

class TestAlbaTraDIS(unittest.TestCase):
	
	def test_small_real(self):
		case = os.path.join(data_dir, 'small_case.insert_site_plot.gz')
		control = os.path.join(data_dir, 'small_control.insert_site_plot.gz')
		emblfile = os.path.join(data_dir, 'annotation.embl')

		t = AlbaTraDIS(TestOptions([case, control], 3, 100, 100, False, 'testoutput', 1, 1, 1, 1, 1, True,1,0, emblfile, 0.1, False, False, 100))
		self.assertTrue(t.run())
		self.assertTrue(os.path.exists('testoutput'))
		self.assertTrue(os.path.exists('testoutput/gene_report.csv'))
		shutil.rmtree("testoutput")
		
	def test_ignore_decreased_insertions(self):
		case = os.path.join(data_dir, 'small_case.insert_site_plot.gz')
		control = os.path.join(data_dir, 'small_control_high_insertions.insert_site_plot.gz')
		emblfile = os.path.join(data_dir, 'annotation.embl')

		t = AlbaTraDIS(TestOptions([case, control], 3, 100, 100, False, 'testoutputx', 1, 1, 1, 1, 1, False,1,0, emblfile, 0.9999, False, False, 100))
		self.assertTrue(t.run())
		self.assertFalse(os.path.exists('.output.pdf'))
		self.assertFalse(os.path.exists('.output.csv'))
		
		self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'expected_no_decrease.plot'), os.path.join('testoutputx', 'combined.plot') ))
		shutil.rmtree("testoutputx")
		
	def test_small_use_annotation(self):
		case = os.path.join(data_dir, 'small_case.insert_site_plot.gz')
		control = os.path.join(data_dir, 'small_control.insert_site_plot.gz')
		emblfile = os.path.join(data_dir, 'annotation.embl')
    
		t = AlbaTraDIS(TestOptions([case, control], 3, 100, 100, False, 'testoutput', 1, 1, 1, 1, 1, True,1,0, emblfile, 0.1, False, True, 100))
		self.assertTrue(t.run())
		self.assertTrue(os.path.exists('testoutput'))
		shutil.rmtree("testoutput")
