'''Given a set of plot files, find the one with the highest number of reads and normalise all the rest into new files'''

from tradistron.PlotParser import PlotParser
from tradistron.PlotGenerator import PlotGenerator
from tempfile import mkstemp

class NormalisePlots:
	
	def __init__(self, plotfiles):
		self.plotfiles = plotfiles
		
	def create_normalised_files(self):
		plot_objs = self.normalise()
		output_files = []
		for p in self.plotfiles:
			fd, output_filename = mkstemp()
			pg = PlotGenerator(plot_objs[p].forward, plot_objs[p].reverse, output_filename)
			pg.construct_file()
			output_files.append(output_filename)
		return output_files
		
	def normalise(self):
		plot_objs = self.read_plots()
		max_plot_reads = self.max_reads(plot_objs)
		
		for p in self.plotfiles:
			current_plot_reads = plot_objs[p].total_reads
			scaling_factor = max_plot_reads/current_plot_reads
			for i,insertion_reads in enumerate(plot_objs[p].forward):
				plot_objs[p].forward[i] = int(plot_objs[p].forward[i]*scaling_factor)
			for i,insertion_reads in enumerate(plot_objs[p].reverse):
				plot_objs[p].reverse[i] = int(plot_objs[p].reverse[i]*scaling_factor)
				
		return plot_objs

	def read_plots(self):
		plot_objs = {}
		for p in self.plotfiles:
			pp = PlotParser(p)
			plot_objs[p] = pp
		return plot_objs
		
	def max_reads(self, plot_objs):
		return max([plot_objs[p].total_reads for p in plot_objs])
