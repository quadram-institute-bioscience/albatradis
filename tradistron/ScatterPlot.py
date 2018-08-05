import numpy
import itertools
import pandas
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
		df = pandas.DataFrame(self.plot_pairs_scatter_coords(self.conditions_plot_objs,1), columns=['rep1','rep2','test_type','coord'], dtype=int)

		ax1 = df.plot.line(x='coord',y='rep1', color= 'red',  label='Condition Rep1')
		df.plot.line(x='coord',y='rep2', color= 'blue',  label='Condition Rep2', ax=ax1)
		
		df = pandas.DataFrame(self.plot_pairs_scatter_coords(self.conditions_plot_objs,2), columns=['rep1','rep2','test_type','coord'], dtype=int)
		df.plot.line(x='coord',y='rep1', color= 'green',  label='Control Rep1', ax=ax1)
		df.plot.line(x='coord',y='rep2', color= 'yellow',  label='Control Rep2', ax=ax1)

		plt.savefig(self.output_filename +"_linear.png" , dpi=100)
		plt.close()
		return self

	
	def create_scatter_plot(self):
		df = pandas.DataFrame(self.plot_pairs_scatter_coords(self.conditions_plot_objs,1), columns=['rep1','rep2','test_type','coord'], dtype=int)
		ax1 = df.plot.scatter(x='rep1', y='rep2',color= 'r',  label='Condition')

		df = pandas.DataFrame(self.plot_pairs_scatter_coords(self.controls_plot_objs,2), columns=['rep1','rep2','test_type','coord'], dtype=int)
		#logx=True, logy=True , loglog=True
		df.plot.scatter(x='rep1', y='rep2',color = 'b',   label='Control', ax=ax1)
		plt.savefig(self.output_filename +"_scatter.png", dpi=100)
		plt.close()
		
		return self
		
	def plot_pairs_scatter_coords(self, plot_objs, type_int ):
		all_coordsout = []
		for plot_obj_pair in itertools.combinations(plot_objs,2):
			window_counts = [self.windows_count(plot_objs[p]) for p in plot_obj_pair] 
			type_values = numpy.full( (1, len(window_counts[0])), type_int)
			range_vals = [numpy.arange(1, len(window_counts[0]) +1, dtype = int)]

			coords = numpy.array(window_counts )
			transformed_coords = numpy.append(coords, type_values, axis=0)
			transformed_coords_index = numpy.append(transformed_coords, range_vals, axis=0).T
			
			if len(all_coordsout) == 0:
				all_coordsout = transformed_coords_index
			else:
				all_coordsout = numpy.append(all_coordsout, transformed_coords_index, axis=0)
		return all_coordsout
		
	def windows_count(self, plot_obj):
		plot_windows = numpy.array_split(plot_obj.combined, self.num_windows)
		return [numpy.sum(p) for p in plot_windows]
