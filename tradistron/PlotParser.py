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
		self.forward = numpy.zeros(self.genome_length )
		self.reverse = numpy.zeros(self.genome_length )
		self.combined = numpy.zeros(self.genome_length )
		
		for i,l in enumerate(lines):	
			insertion_count_clean = l.split()

			f = int(float(insertion_count_clean[0]))
			r = int(float(insertion_count_clean[1]))
			if numpy.absolute(f) >= self.minimum_threshold:
				self.forward[i] = f
				
			if numpy.absolute(r) >= self.minimum_threshold:
				self.reverse[i] = r

			if numpy.absolute(f) + numpy.absolute(r) >= self.minimum_threshold:
				self.combined[i] = f + r

		return self
			
