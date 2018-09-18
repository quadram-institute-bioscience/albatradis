'''Given a set of plot files, find the one with the highest number of reads and normalise all the rest into new files'''

from albatradis.PlotParser import PlotParser
from albatradis.PlotGenerator import PlotGenerator
from tempfile import mkstemp
import numpy

class NormalisePlots:
	
	def __init__(self, plotfiles, minimum_proportion_insertions, verbose = False):
		self.plotfiles = plotfiles
		self.verbose = verbose
		self.minimum_proportion_insertions = minimum_proportion_insertions
		self.plot_objs = self.read_plots()
		
	def create_normalised_files(self):
		plot_objs = self.normalise()
		output_files = []
		for p in self.plotfiles:
			fd, output_filename = mkstemp()
			pg = PlotGenerator(plot_objs[p].forward, plot_objs[p].reverse, output_filename)
			pg.construct_file()
			output_files.append(output_filename)
		return output_files
		
	def decreased_insertion_reporting(self):
		max_plot_reads = self.max_reads(self.plot_objs)
		for t in self.plot_total_reads(self.plot_objs):
			if t/max_plot_reads < self.minimum_proportion_insertions:
				print("No. of reads in file is "+str(t)+ " compared to a maximum of "+str(max_plot_reads) + " so we cant call decreased insertions accurately")
				return False
				
		# check number of insertion sites
		max_plot_insertions = self.max_insertions (self.plot_objs)
		for t in self.plot_insertions(self.plot_objs):
			if t/max_plot_insertions < self.minimum_proportion_insertions:
				print("No. of insertions in file is "+str(t)+ " compared to a maximum of "+str(max_plot_insertions) + " so we cant call decreased insertions accurately")
				return False
		return True
		
	def normalise(self):
		max_plot_reads = self.max_reads(self.plot_objs)
		
		for p in self.plotfiles:
			current_plot_reads = self.plot_objs[p].total_reads
			scaling_factor = max_plot_reads/current_plot_reads
			if self.verbose:
				print("\t".join(("Normalise", p, str(current_plot_reads), str(max_plot_reads), str(scaling_factor))))
			
			self.plot_objs[p].forward = numpy.multiply(self.plot_objs[p].forward, scaling_factor, dtype = float, casting = 'unsafe')
			self.plot_objs[p].reverse = numpy.multiply(self.plot_objs[p].reverse, scaling_factor, dtype = float, casting = 'unsafe')
				
		return self.plot_objs

	def read_plots(self):
		plot_objs = {}
		for p in self.plotfiles:
			pp = PlotParser(p)
			plot_objs[p] = pp
		return plot_objs
		
	def plot_insertions(self, plot_objs):
		ins = [plot_objs[p].total_insertions for p in plot_objs]
		return ins
		
	def plot_total_reads(self, plot_objs):
		reads = [plot_objs[p].total_reads for p in plot_objs]
		return reads
	
	def max_reads(self, plot_objs):
		return max(self.plot_total_reads(plot_objs))
		
	def max_insertions(self, plot_objs):
		return max(self.plot_insertions(plot_objs))
