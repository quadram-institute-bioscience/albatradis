'''Take in the input files, parse them to create new files.'''
from tempfile import mkstemp
from albatradis.EMBLGenerator import EMBLGenerator
from albatradis.EMBLExpandGenes import EMBLExpandGenes
from albatradis.PlotParser import PlotParser
from albatradis.WindowGenerator import WindowGenerator
from albatradis.PlotGenerator import PlotGenerator

class PrepareEMBLFile:
	def __init__(self, plotfile, minimum_threshold, window_size, window_interval, use_annotation, prime_feature_size,emblfile):
		self.plotfile = plotfile
		self.minimum_threshold = minimum_threshold
		self.window_size = window_size
		self.window_interval = window_interval
		self.use_annotation    = use_annotation
		self.prime_feature_size = prime_feature_size
		self.emblfile = emblfile
		
		self.forward_plot_filename = ""
		self.reverse_plot_filename = ""
		self.combined_plot_filename = ""
		self.embl_filename = ""
		
	def genome_length(self):
		return self.plot_parser_obj.genome_length
		
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

	def embl_file_expand_genes(self):
		fd, embl_filename = mkstemp()
		
		eg = EMBLExpandGenes(self.emblfile, self.prime_feature_size)
		eg.construct_file(embl_filename)
		return embl_filename
		
	def create_file(self):
		self.plot_parser_obj = self.plot_parser()
		
		# use the annotation and add 3/5 prime blocks to each gene
		if self.use_annotation:
			self.embl_filename = self.embl_file_expand_genes()
		else:			
			# Use sliding windows only
			self.embl_filename = self.create_embl_file()

		return self.embl_filename		
