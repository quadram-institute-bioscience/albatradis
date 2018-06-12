'''Takes in arrays for forward and reverse integers and creates a new file'''
class PlotGenerator:
	def __init__(self, forward,reverse, filename):
		self.forward = forward
		self.reverse = reverse
		self.filename = filename
		
		
	def plot_length(self):
		total_length = len(self.forward)
		if len(self.reverse) > total_length:
			total_length = len(self.reverse) 
			
		return total_length
		
	def construct_line(self, i):
		line = ""
		if len(self.forward) > i:
			line += str(self.forward[i])+ " "
		else:
			line += "0 "
		
		if len(self.reverse) > i:
			line += str(self.reverse[i])+ "\n"
		else:
			line += "0\n"
		return line
		
	def construct_file(self):
		with open(self.filename, 'w') as plotfile:
			for i in range(0,self.plot_length()):
				plotfile.write(self.construct_line(i))
		return self