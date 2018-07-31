import os
import sys
import re
import numpy
import subprocess
import pandas

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

class PlotParser:
	
	def __init__(self, filename, minimum_threshold = 0):
		self.filename = filename
		self.minimum_threshold = minimum_threshold
		self.forward = []
		self.reverse = []
		self.combined = []
		
		self.split_lines()
		
		self.total_reads = sum(self.combined)
		self.total_insertions = sum([1 for a in self.combined if a > 0 ])
		
	def split_lines(self):
		insert_site_array = pandas.read_csv(self.filename, delim_whitespace=True, dtype=int, engine='c', header=None).values

		self.genome_length =  len(insert_site_array)
		
		self.forward = insert_site_array[:,0]
		self.reverse = insert_site_array[:,1]
		
		if self.minimum_threshold != 0:
			self.forward = self.filter_column(self.forward, self.genome_length)
			self.reverse = self.filter_column(self.reverse, self.genome_length)
		
		self.combined = [ self.forward[i] + self.reverse[i] for i in range(0, self.genome_length)]
		
		return self
	
	def filter_column(self,ins_array,genome_length ):
		abs_ins_values = numpy.absolute(ins_array)
		return [ ins_array[i] if  abs_ins_values[i] >= self.minimum_threshold else 0  for i in range(0, genome_length)]
		
