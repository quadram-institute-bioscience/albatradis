'''Driver class'''
import logging
import os
import sys

from tradistron.ScatterPlot import ScatterPlot

class TradisTronScatterPlot:
	def __init__(self, options):
		self.logger            = logging.getLogger(__name__)
		self.window_size       = options.window_size
		self.verbose           = options.verbose
		self.outputfile        = options.outputfile
		self.controls          = options.control
		self.conditions        = options.conditions

		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)


	def run(self):
		p = ScatterPlot(self.conditions, self.controls, self.window_size, self.outputfile, verbose = self.verbose)
		p.create_scatter_plot()
		return self
