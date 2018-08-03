import numpy
import itertools
import pandas
import matplotlib.pyplot as plt
from tradistron.PlotParser import PlotParser

class ScatterPlot:
	#assumption is that there are 2 or more conditions, and 2 or more controls.
	def __init__(self, conditions, controls, window_size, output_filename):
		self.conditions = conditions
		self.controls = controls
		self.window_size = window_size
		self.output_filename = output_filename
		
		self.conditions_plot_objs = self.get_plot_objs(self.conditions)
		self.controls_plot_objs = self.get_plot_objs(self.controls)
		
	def get_plot_objs(self, files):
		plot_objs = {}
		for f in files:
			plot_objs[f] = PlotParser(f)
		return plot_objs
		
	def create_scatter_plot(self):
		# 0 = conditions
		# 1 = controls
		coords = self.plot_pairs_scatter_coords(self.conditions_plot_objs,0) + self.plot_pairs_scatter_coords(self.controls_plot_objs,1)
		df = pandas.DataFrame(coords, columns=['rep1','rep2','type'])
		df.plot.scatter(x='rep1', y='rep2', c='type')
		df.plot()
		plt.show()
		
		return self
		
	def plot_pairs_scatter_coords(self, plot_objs, type_int ):
		all_coordsout = []
		for plot_obj_pair in itertools.combinations(plot_objs,2):
			window_counts = [self.windows_count(plot_objs[p]) for p in plot_obj_pair] 
			type_values = numpy.full( (1, len(window_counts[0])), type_int)
			coords = numpy.array(window_counts )
			transformed_coords = numpy.append(coords, type_values, axis=0).T
			print(transformed_coords)
			
			if len(all_coordsout) == 0:
				all_coordsout = transformed_coords
			else:
				all_coordsout = numpy.append(all_coordsout, transformed_coords, axis=0)
		return all_coordsout
		
	def windows_count(self, plot_obj):
		plot_windows = numpy.array_split(plot_obj.combined, self.window_size)
		return [numpy.sum(p) for p in plot_windows]
