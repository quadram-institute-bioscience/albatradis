class Window:
	def __init__(self, start, end):
		self.start = start
		self.end = end
		
	def name(self):
		return str(self.start) + "_" + str(self.end)