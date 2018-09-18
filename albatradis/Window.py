class Window:
	def __init__(self, start, end):
		self.start = start
		self.end = end
		
	def name_start_one(self):
		return str(self.start + 1) + "_" + str(self.end)
