from albatradis.PlotParser import PlotParser
from albatradis.PlotGenerator import PlotGenerator
from tempfile import mkstemp

class PlotMasking:
	def __init__(self, insertion_plot_files, masking_plot_file):
		self.insertion_plot_files   = insertion_plot_files
		self.masking_plot_file 		= masking_plot_file
		self.output_plot_files      = self.mask_plot_files()

	'''Take in a logfc plot file and apply it as a mask to the insertion plot files'''
	def mask_plot_files(self):
		masking_plot = PlotParser(self.masking_plot_file)
		output_plot_files = {}
		
		for insertion_filename in self.insertion_plot_files:
			insertion_plot = PlotParser(insertion_filename)
			
			for i, mask in enumerate(masking_plot.combined):
				if mask == 0:
					insertion_plot.forward[i] = 0
					insertion_plot.reverse[i] = 0
					
			fd, output_filename = mkstemp()
			pg = PlotGenerator(insertion_plot.forward, insertion_plot.reverse, output_filename)
			pg.construct_file()
			output_plot_files[insertion_filename] = output_filename
		return output_plot_files
