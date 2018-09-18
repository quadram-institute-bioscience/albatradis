import numpy
import itertools
import pandas
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from tradistron.PlotParser import PlotParser
from tradistron.NormalisePlots import NormalisePlots

class ScatterPlot:
	#assumption is that there are 2 or more conditions, and 2 or more controls.
	def __init__(self, conditions, controls, window_size, output_filename,normalise,  verbose = False):
		self.conditions = conditions
		self.controls = controls
		self.window_size = window_size
		self.output_filename = output_filename
		self.verbose = verbose
		self.normalise = normalise
		
		if normalise:
			n = NormalisePlots(self.conditions + self.controls, 0.0000001, verbose = self.verbose)
			plotfiles = n.create_normalised_files()
			self.conditions = plotfiles[0:len(self.conditions)]
			self.controls = plotfiles[len(self.conditions):]
		
		self.conditions_plot_objs = self.get_plot_objs(self.conditions)
		self.controls_plot_objs = self.get_plot_objs(self.controls)
		self.set_num_windows()
		
	
	def set_num_windows(self):
			self.num_windows = numpy.ceil(self.genome_size/self.window_size)
	
	def get_plot_objs(self, files):
		plot_objs = {}
		for f in files:
			plot_objs[f] = PlotParser(f)
			self.genome_size = len(plot_objs[f].combined)
			if self.verbose:
				print(plot_objs[f])
		return plot_objs
		
	def create_linear_plot(self):
		cond = self.plot_pairs_scatter_coords(self.conditions_plot_objs)
		cont = self.plot_pairs_scatter_coords(self.controls_plot_objs)

		df = pandas.DataFrame({'Condition-Rep1': cond[:,0], 'Condition-Rep2': cond[:,1], 'Control-Rep1': cont[:,0], 'Control-Rep2': cont[:,1]}, index = cond[:,2].astype(int))
		
		df = df.astype(int)

		ax1 = df.plot.line(colormap='summer')
		ax1.set_xlabel("Genome (base position)")
		ax1.set_ylabel("No. of insertions")
		plt.title("Insertions binned into Windows of "+str(self.window_size)+" bases")
		plt.savefig(self.output_filename +"_linear.png" , dpi=100)
		
		plt.close()
		return self
	
	def create_scatter_plot(self):
		cond = self.plot_pairs_scatter_coords(self.conditions_plot_objs)
		cont = self.plot_pairs_scatter_coords(self.controls_plot_objs)
		
		df = pandas.DataFrame({'Condition-Rep1': cond[:,0], 'Condition-Rep2': cond[:,1], 'Control-Rep1': cont[:,0], 'Control-Rep2': cont[:,1]}, index = cond[:,1].astype(int))
		df = df.astype(int)
		
		ax1 = df.plot.scatter(x='Condition-Rep1', y='Condition-Rep2',color= 'r',  label='Condition',  loglog=True)
		df.plot.scatter(x='Control-Rep1', y='Control-Rep2',color = 'b',   label='Control', ax=ax1, loglog=True)
		
		ax1.set_xlabel("No. of insertions (Rep 2)")
		ax1.set_ylabel("No. of insertions (Rep 1)")
		plt.title("No of insertions in Rep 1 vs Rep 2 for a window of "+str(self.window_size)+" bases")
		plt.savefig(self.output_filename +"_scatter.png", dpi=100)
		plt.close()
		
		return self
		
		
	def create_abs_scatter_plot(self):
		cond = self.abs_change_axis(self.plot_pairs_scatter_coords(self.conditions_plot_objs))
		cont = self.abs_change_axis(self.plot_pairs_scatter_coords(self.controls_plot_objs))
		
		df = pandas.DataFrame({'Condition-Rep1': cond[:,0], 'Condition-Rep2': cond[:,1], 'Control-Rep1': cont[:,0], 'Control-Rep2': cont[:,1]}, index = cond[:,1].astype(int))
		df = df.astype(int)
		
		ax1 = df.plot.scatter(x='Condition-Rep2', y='Condition-Rep1',color= 'r',  label='Condition',  loglog=True)
		df.plot.scatter(x='Control-Rep2', y='Control-Rep1',color = 'b',   label='Control', ax=ax1, loglog=True)
		
		ax1.set_xlabel("No. of insertions (Rep 2)")
		ax1.set_ylabel("Normalised No. of insertions (Rep 1)")
		plt.title("Normalised insertions Rep 1 vs Rep 2 for a window of "+str(self.window_size)+" bases")
		
		plt.savefig(self.output_filename +"_absscatter.png", dpi=100)
		plt.close()
		return self
		
	def abs_change_axis(self, coords):
		abs_vals = numpy.array([numpy.absolute(numpy.subtract(coords[:,1],coords[:,0]), dtype = int)])
		index_array = numpy.array([coords[:,0]])
		return numpy.append(abs_vals,index_array, axis=0).T

		
	def plot_pairs_scatter_coords(self, plot_objs ):
		all_coordsout = []
		for plot_obj_pair in itertools.combinations(plot_objs,2):
			window_counts = [self.windows_count(plot_objs[p]) for p in plot_obj_pair] 
			range_vals = [numpy.multiply(self.window_size, numpy.arange(0, len(window_counts[0]) ), dtype = int)]

			coords = numpy.array(window_counts )
			transformed_coords_index = numpy.append(coords, range_vals, axis=0).T
			
			if len(all_coordsout) == 0:
				all_coordsout = transformed_coords_index
			else:
				all_coordsout = numpy.append(all_coordsout, transformed_coords_index, axis=0)
		return all_coordsout
		
	def windows_count(self, plot_obj):
		plot_windows = numpy.array_split(plot_obj.combined, self.num_windows)
		return [numpy.sum(p) for p in plot_windows]
