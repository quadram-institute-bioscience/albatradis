'''Driver class'''
import logging
import os
import sys

from albatradis.ScatterPlot import ScatterPlot

class AlbaTraDISScatterPlot:
	def __init__(self, options):
		self.logger            = logging.getLogger(__name__)
		self.window_size       = options.window_size
		self.verbose           = options.verbose
		self.outputfile        = options.outputfile
		self.controls          = options.control
		self.conditions        = options.condition
		self.normalise         = options.normalise

		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)


	def run(self):
		p = ScatterPlot(self.conditions, self.controls, self.window_size, self.outputfile, self.normalise, verbose = self.verbose)
		
		# High level windows
		for ws in [10000,100000]:
			if self.verbose:
				print("Running Sliding window size:\t"+str(ws))
			p.output_filename = self.outputfile+'_'+str(ws)
			p.window_size = ws
			p.set_num_windows()
			p.create_scatter_plot()
			p.create_linear_plot()
			p.create_abs_scatter_plot()

		p.output_filename = self.outputfile
		p.window_size = self.window_size
		p.set_num_windows()
		p.create_scatter_plot()
		p.create_linear_plot()
		p.create_abs_scatter_plot()
		return self
