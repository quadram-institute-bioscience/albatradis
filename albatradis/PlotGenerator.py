import pyximport
pyximport.install()

from file_manipulation import write_plot_file

import numpy
import pandas

'''Takes in arrays for forward and reverse integers and creates a new file'''
class PlotGenerator:
	def __init__(self, forward,reverse, filename):
		self.forward = forward
		self.reverse = reverse
		self.filename = filename
		
		self.forward_length = len(self.forward)
		self.reverse_length = len(self.reverse)
		
	def construct_file(self):

		p_len = self.plot_length()
		if len(self.forward) == 0:
			self.forward = numpy.zeros( p_len, dtype=float )

		if len(self.reverse) == 0:
			self.reverse = numpy.zeros( p_len, dtype=float )
			
		write_plot_file(self.filename, self.forward, self.reverse, p_len)	

		return self
		
	def plot_length(self):
		total_length = self.forward_length
		if self.reverse_length > total_length:
			total_length = self.reverse_length
			
		return total_length
