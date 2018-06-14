import csv
import numpy
from tempfile import mkstemp

class LogFC:
	def __init__(self, start, end, logfc_value):
		self.start = start
		self.end = end
		self.logfc_value = logfc_value

class PlotLog:
	def __init__(self, comparison_filename, genome_length, minimum_logfc, pvalue, minimum_logcpm):
		self.comparison_filename = comparison_filename
		self.genome_length = genome_length
		self.minimum_logfc = minimum_logfc
		self.pvalue = pvalue
		self.minimum_logcpm = minimum_logcpm
		
		fd, self.output_filename = mkstemp()

	def construct_plot_file(self):
		logfc_coord_values = self.read_comparison_file()
		logfc_to_bases = self.genome_wide_logfc(logfc_coord_values)
		
		with open(self.output_filename, 'w') as plotfile:
			for i in logfc_to_bases:
				plotfile.write(self.construct_line(i))
		return self
		
	def construct_line(self, i):
		line = ""
		if i >= 0:
			line += str(int(i))+ "\t"
		else:
			line += "0\t"
		
		if i < 0:
			line += str(int(i))+ "\n"
		else:
			line += "0\n"

		return line
		
	def genome_wide_logfc(self,logfc_coord_values):
		logfc_to_bases = numpy.zeros(self.genome_length)
		
		for l in logfc_coord_values:
			for i in range(l.start -1, l.end):
				
				if logfc_to_bases[i] < numpy.absolute(l.logfc_value):
					logfc_to_bases[i] = l.logfc_value
		return logfc_to_bases
			
		
	def read_comparison_file(self):
		logfc_coord_values = []
		
		with open(self.comparison_filename, newline='') as csvfile:
			comparison_reader = csv.reader(csvfile, delimiter=',')
			for row in comparison_reader:
				logfc = row[3]
				if logfc == 'logFC':
					 continue
					 
				logfc = int(float(row[3]))
				
				if numpy.absolute(logfc) < self.minimum_logfc or int(float(row[5])) >= self.pvalue:
					logfc = 0
					
				if numpy.absolute(int(float(row[4]))) < self.minimum_logcpm:
					logfc = 0
					
				# encodes coordinates
				gene_name = row[0]
				start, end = gene_name.split("_")
				
				logfc_coord_values.append(LogFC(int(start), int(end), logfc))
				
		return logfc_coord_values