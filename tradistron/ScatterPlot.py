import numpy
import itertools
import pandas
import matplotlib.pyplot as plt
from tradistron.PlotParser import PlotParser

class ScatterPlot:
	#assumption is that there are 2 or more conditions, and 2 or more controls.
	def __init__(self, conditions, controls, window_size, output_filename, verbose = False):
		self.conditions = conditions
		self.controls = controls
		self.window_size = window_size
		self.output_filename = output_filename
		self.verbose = verbose
		
		self.conditions_plot_objs = self.get_plot_objs(self.conditions)
		self.controls_plot_objs = self.get_plot_objs(self.controls)
		self.num_windows = numpy.ceil(self.genome_size/self.window_size)
		
	def get_plot_objs(self, files):
		plot_objs = {}
		for f in files:
			plot_objs[f] = PlotParser(f)
			self.genome_size = len(plot_objs[f].combined)
			if self.verbose:
				print(plot_objs[f])
		return plot_objs
		
	def create_scatter_plot(self):
		# 0 = conditions
		# 1 = controls
		df = pandas.DataFrame(self.plot_pairs_scatter_coords(self.conditions_plot_objs,1), columns=['rep1','rep2','test_type'])
		ax1 = df.plot.scatter(x='rep1', y='rep2',color= 'r',  label='Condition')

		df = pandas.DataFrame(self.plot_pairs_scatter_coords(self.controls_plot_objs,2), columns=['rep1','rep2','test_type'])
		#logx=True, logy=True , loglog=True
		df.plot.scatter(x='rep1', y='rep2',color = 'b',   label='Control', ax=ax1)

		plt.savefig(self.output_filename , dpi=100)
		plt.close()
		
		return self
		
	def plot_pairs_scatter_coords(self, plot_objs, type_int ):
		all_coordsout = []
		for plot_obj_pair in itertools.combinations(plot_objs,2):
			window_counts = [self.windows_count(plot_objs[p]) for p in plot_obj_pair] 
			type_values = numpy.full( (1, len(window_counts[0])), type_int)
			coords = numpy.array(window_counts )
			transformed_coords = numpy.append(coords, type_values, axis=0).T
			
			if len(all_coordsout) == 0:
				all_coordsout = transformed_coords
			else:
				all_coordsout = numpy.append(all_coordsout, transformed_coords, axis=0)
		return all_coordsout
		
	def windows_count(self, plot_obj):
		plot_windows = numpy.array_split(plot_obj.combined, self.num_windows)
		return [numpy.sum(p) for p in plot_windows]
