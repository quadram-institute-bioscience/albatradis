'''Take in the input files, parse them to create new files.'''
from tempfile import mkstemp
from tradistron.EMBLGenerator import EMBLGenerator
from tradistron.PlotParser import PlotParser
from tradistron.WindowGenerator import WindowGenerator
from tradistron.PlotGenerator import PlotGenerator

class PrepareInputFiles:
	def __init__(self, plotfile, minimum_threshold, window_size, window_interval):
		self.plotfile = plotfile
		self.minimum_threshold = minimum_threshold
		self.window_size = window_size
		self.window_interval = window_interval
		
		self.forward_plot_filename = ""
		self.reverse_plot_filename = ""
		self.combined_plot_filename = ""
		self.embl_filename = ""
		
	def plot_parser(self):
		return PlotParser(self.plotfile,self.minimum_threshold)
		
	def windows(self):
		w = WindowGenerator(self.plot_parser_obj.genome_length, self.window_size, self.window_interval)
		return w.create_windows()
		
	def create_embl_file(self):
		e = EMBLGenerator(self.windows(), self.plot_parser_obj.genome_length)

		fd, embl_filename = mkstemp()
		e.construct_file(embl_filename)
		return embl_filename

	def create_split_plot_file(self, forward, reverse):
		fd, filename = mkstemp()
		p = PlotGenerator(forward, reverse, filename)
		p.construct_file()
		return filename
		
	def create_all_files(self):
		self.plot_parser_obj = self.plot_parser()
		self.embl_filename = self.create_embl_file()
		
		self.forward_plot_filename = self.create_split_plot_file(self.plot_parser_obj.forward, [])
		self.reverse_plot_filename = self.create_split_plot_file([], self.plot_parser_obj.reverse)
		self.combined_plot_filename = self.create_split_plot_file(self.plot_parser_obj.forward,  self.plot_parser_obj.reverse)
		return self
		
		
