'''Take in the input files, parse them to create new files.'''
from tempfile import mkstemp
from albatradis.EMBLGenerator import EMBLGenerator
from albatradis.EMBLExpandGenes import EMBLExpandGenes
from albatradis.PlotParser import PlotParser
from albatradis.WindowGenerator import WindowGenerator
from albatradis.PlotGenerator import PlotGenerator

class PrepareInputFiles:
	def __init__(self, plotfile, minimum_threshold):
		self.plotfile = plotfile
		self.minimum_threshold = minimum_threshold
		
		self.forward_plot_filename = ""
		self.reverse_plot_filename = ""
		self.combined_plot_filename = ""
		self.embl_filename = ""
		
	def genome_length(self):
		return self.plot_parser_obj.genome_length
		
	def plot_parser(self):
		return PlotParser(self.plotfile,self.minimum_threshold)

	def create_split_plot_file(self, forward, reverse):
		fd, filename = mkstemp()
		p = PlotGenerator(forward, reverse, filename)
		p.construct_file()
		return filename
		
	def create_all_files(self):
		self.plot_parser_obj = self.plot_parser()
		
		self.forward_plot_filename = self.create_split_plot_file(self.plot_parser_obj.forward, [])
		self.reverse_plot_filename = self.create_split_plot_file([], self.plot_parser_obj.reverse)
		self.combined_plot_filename = self.create_split_plot_file(self.plot_parser_obj.forward,  self.plot_parser_obj.reverse)
		return self		
