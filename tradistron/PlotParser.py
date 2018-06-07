import os
import sys
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
		
		self.genome_length = len(self.combined)
		

	def get_filehandle(self):
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

		return lines
		
	def split_lines(self):
		for i,l in enumerate(self.read_file()):
			insertion_count = l.split()
			
			# sanitise input
			insertion_count_clean = [s for s in insertion_count if s.isdigit()]
			if len(insertion_count_clean) != 2:
				raise InvalidFileFormat("Invalid line in file: " + str(l))
			
			if int(insertion_count_clean[0]) >= self.minimum_threshold:
				self.forward.append(int(insertion_count_clean[0])) 
			else:
				self.forward.append(0)
				
			if int(insertion_count_clean[1]) >= self.minimum_threshold:
				self.reverse.append(int(insertion_count_clean[1])) 
			else:
				self.reverse.append(0)
			
			if int(insertion_count_clean[0]) + int( insertion_count_clean[1]) >= self.minimum_threshold:
				self.combined.append(int(insertion_count_clean[0]) + int( insertion_count_clean[1]))
			else:
				self.combined.append(0)
		return self
			
