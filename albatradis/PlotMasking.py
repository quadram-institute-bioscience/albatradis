from albatradis.PlotParser import PlotParser
from albatradis.PlotGenerator import PlotGenerator
from tempfile import mkstemp
import shutil

class PlotMasking:
	def __init__(self, insertion_plot_files, masking_plot_file, strict_signal):
		self.insertion_plot_files   = insertion_plot_files
		self.masking_plot_file 		= masking_plot_file
		self.strict_signal          = strict_signal
		self.output_plot_files      = self.mask_plot_files()

	'''Take in a logfc plot file and apply it as a mask to the insertion plot files'''
	def mask_plot_files(self):
		if self.strict_signal:
			return self.masking()
		else:
			return self.no_masking()
		
	def masking(self):
		masking_plot = PlotParser(self.masking_plot_file)
		output_plot_files = {}
		
		for insertion_filename in self.insertion_plot_files:
			insertion_plot = PlotParser(insertion_filename)
			
			for i, mask in enumerate(masking_plot.combined):
				if mask == 0 and self.strict_signal:
					insertion_plot.forward[i] = 0
					insertion_plot.reverse[i] = 0
					
			fd, output_filename = mkstemp()
			pg = PlotGenerator(insertion_plot.forward, insertion_plot.reverse, output_filename)
			pg.construct_file()
			output_plot_files[insertion_filename] = output_filename
		return output_plot_files
		
	def no_masking(self):
		output_plot_files = {}
		for insertion_filename in self.insertion_plot_files:
			fd, output_filename = mkstemp()
			shutil.copy(insertion_filename,output_filename)
			output_plot_files[insertion_filename] = output_filename
		return output_plot_files
