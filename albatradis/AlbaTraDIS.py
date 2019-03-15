'''Driver class'''
import logging
import os
import sys
import time
import re

from albatradis.BlockInsertions import BlockInsertions
from albatradis.NormalisePlots import NormalisePlots


class AlbaTraDIS:
	def __init__(self, options):
		self.logger            = logging.getLogger(__name__)
		self.plotfiles         = options.plotfiles
		self.minimum_threshold = options.minimum_threshold
		self.window_size       = options.window_size
		self.window_interval   = options.window_interval
		self.verbose           = options.verbose
		self.minimum_logfc     = options.minimum_logfc
		self.pvalue            = options.pvalue
		self.prefix            = options.prefix
		self.minimum_logcpm    = options.minimum_logcpm
		self.iterations        = options.iterations
		self.dont_normalise_plots   = options.dont_normalise_plots
		self.minimum_block     = options.minimum_block
		self.span_gaps         = options.span_gaps
		self.emblfile          = options.emblfile
		self.minimum_proportion_insertions = options.minimum_proportion_insertions
		self.strict_signal     = options.strict_signal
		self.use_annotation    = options.use_annotation
		self.prime_feature_size = options.prime_feature_size
		
		self.genome_length = 0
		
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
			
		self.blocks = []
	
	def run(self):
		plotfiles = self.plotfiles
		report_decreased_insertions = True
		if not self.dont_normalise_plots:
			n = NormalisePlots(self.plotfiles, self.minimum_proportion_insertions)
			plotfiles = n.create_normalised_files()
			report_decreased_insertions = n.decreased_insertion_reporting()
		
		if self.iterations == 1:
			bi = BlockInsertions(self.logger, plotfiles, self.minimum_threshold, self.window_size, self.window_interval, self.verbose, self.minimum_logfc, self.pvalue, self.prefix, self.minimum_logcpm, self.minimum_block, self.span_gaps, self.emblfile, report_decreased_insertions,self.strict_signal,self.use_annotation, self.prime_feature_size)
			bi.run()
			self.blocks = bi.blocks
			plotfiles = bi.output_plots.values()
		else:
		
			for i in range(1,self.iterations+1):
				bi = BlockInsertions(self.logger, plotfiles, self.minimum_threshold, self.window_size, self.window_interval, self.verbose, self.minimum_logfc, self.pvalue, self.prefix + "_" +str(i), self.minimum_logcpm, self.minimum_block, self.span_gaps, self.emblfile, report_decreased_insertions,self.strict_signal,self.use_annotation, self.prime_feature_size)
				bi.run()
				self.blocks = bi.blocks
				plotfiles = bi.output_plots.values()
		
		self.cleanup()
		return self
		
	def cleanup(self):
		for f in ['.output.pdf', '.output.csv']:
			if os.path.exists(f):
				os.remove(f)
				
		# remove any tmp files left over
		for dirname, dirnames, filenames in os.walk(self.prefix):
			for filename in filenames:
				deletion_regex = "^tmp........$"
				full_path_of_file_for_deletion = os.path.join(self.prefix, filename)
				if(re.match(str(deletion_regex), filename) != None and os.path.exists(full_path_of_file_for_deletion)):
					if self.verbose > 0:
						print("Deleting file: "+ full_path_of_file_for_deletion + " regex:"+deletion_regex)
					os.remove(full_path_of_file_for_deletion)
		
