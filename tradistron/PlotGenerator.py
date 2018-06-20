import numpy

'''Takes in arrays for forward and reverse integers and creates a new file'''
class PlotGenerator:
	def __init__(self, forward,reverse, filename):
		self.forward = forward
		self.reverse = reverse
		self.filename = filename
		
		self.forward_length = len(self.forward)
		self.reverse_length = len(self.reverse)
		
	def create_body(self):
		self.normalised_forward = numpy.zeros( self.plot_length(), dtype=int )
		self.normalised_reverse = numpy.zeros( self.plot_length(), dtype=int )
		
		for i,val in enumerate(self.forward):
			if val != 0:
				self.normalised_forward[i]  = val
				
		for i,val in enumerate(self.reverse):
			if val != 0:
				self.normalised_reverse[i]  = val
		
		return "\n".join([str(self.normalised_forward[i]) + " " + str(self.normalised_reverse[i]) for i in range(self.plot_length()) ]) + "\n"
		
	def plot_length(self):
		total_length = self.forward_length
		if self.reverse_length > total_length:
			total_length = self.reverse_length
			
		return total_length
		
	def construct_file(self):
		with open(self.filename, 'w') as plotfile:
			plotfile.write(self.create_body())
		return self