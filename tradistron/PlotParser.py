import os
import sys
import re
import numpy
import subprocess

class ErrorReadingFile (Exception): pass
class InvalidFileFormat (Exception): pass

class PlotParser:
	
	def __init__(self, filename, minimum_threshold = 1):
		self.filename = filename
		self.minimum_threshold = minimum_threshold
		self.forward = []
		self.reverse = []
		self.combined = []
		
		self.split_lines()
		
		self.total_reads = sum(self.combined)
		self.total_insertions = sum([1 for a in self.combined if a > 0 ])
		

	def get_filehandle(self):
		if not os.path.exists(self.filename):
			raise ErrorReadingFile("Error opening for reading file '" + self.filename + "'")
		
		if self.filename == '-':
			f = sys.stdin
		elif self.filename.endswith('.gz'):
			retcode = subprocess.call('gunzip -t ' + self.filename, shell=True)
			f = os.popen('gunzip -c ' + self.filename)
		else:
			try:
				f = open(self.filename)
			except:
				raise ErrorReadingFile("Error opening for reading file '" + self.filename + "'")
		
		return f

	def read_file(self):
		f = self.get_filehandle()
		lines = [line.strip() for line in f]
		f.close()

		return lines
		
	def split_lines(self):
		lines = self.read_file()
		self.genome_length =  len(lines)
		self.forward = numpy.zeros(self.genome_length, dtype=int )
		self.reverse = numpy.zeros(self.genome_length, dtype=int )
		self.combined = numpy.zeros(self.genome_length, dtype=int )
		
		split_lines = [l.split() for l in lines]
		int_split_lines = [ [int(float(l[0])), int(float(l[1]))] for l in split_lines]

		for i,l in enumerate(int_split_lines):	
			if numpy.absolute(l[0]) >= self.minimum_threshold:
				self.forward[i] = l[0]
				
			if numpy.absolute(l[1]) >= self.minimum_threshold:
				self.reverse[i] = l[1]

			if numpy.absolute(l[0]) + numpy.absolute(l[1]) >= self.minimum_threshold:
				self.combined[i] = l[0] + l[1]

		return self
			
