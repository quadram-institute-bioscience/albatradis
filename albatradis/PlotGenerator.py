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
			self.forward = numpy.zeros( p_len, dtype=int )

		if len(self.reverse) == 0:
			self.reverse = numpy.zeros( p_len, dtype=int )

		df = pandas.DataFrame({'forward': self.forward, 'reverse': self.reverse}, dtype = int)
		
		df.to_csv(self.filename, header=None, index = False, sep = ' ', chunksize=100000)
		return self
		
	def plot_length(self):
		total_length = self.forward_length
		if self.reverse_length > total_length:
			total_length = self.reverse_length
			
		return total_length
