from tradistron.Window import Window

class WindowGenerator:
	def __init__(self, genome_length, window_size, window_interval):
		self.genome_length = genome_length
		self.window_size = window_size
		self.window_interval = window_interval
		
	def create_windows(self):
		if self.genome_length < self.window_size or self.genome_length < self.window_interval or self.window_interval < 1 or self.genome_length  < 1 or self.window_size < 1:
			return []
		
		end = self.genome_length - self.window_size +1
		windows = [ Window(i,i+self.window_size) for i in range(0,end, self.window_interval)]
		return windows
