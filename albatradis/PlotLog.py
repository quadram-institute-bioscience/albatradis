import csv
import numpy
import pandas
from tempfile import mkstemp
from albatradis.Block import Block
from albatradis.EMBLReader import EMBLReader


class LogFC:
	def __init__(self, start, end, logfc_value):
		self.start = start
		self.end = end
		self.logfc_value = logfc_value

class PlotLog:
	def __init__(self, comparison_filename, genome_length, minimum_logfc, pvalue, qvalue, minimum_logcpm, window_size, span_gaps, report_decreased_insertions, embl_file):
		self.comparison_filename = comparison_filename
		self.genome_length  = genome_length
		self.minimum_logfc  = minimum_logfc
		self.pvalue         = pvalue
		self.qvalue         = qvalue
		self.minimum_logcpm = minimum_logcpm
		self.window_size    = window_size
		self.span_gaps      = span_gaps
		self.report_decreased_insertions = report_decreased_insertions
		self.embl_file      = embl_file
		
		fd, self.output_filename = mkstemp()

	def construct_plot_file(self):
		logfc_coord_values = self.read_comparison_file()
		logfc_to_bases = self.genome_wide_logfc(logfc_coord_values)

		forward_logfc = [i if i >= 0 else 0.0 for i in logfc_to_bases]
		reverse_logfc = [i if i < 0 else 0.0 for i in logfc_to_bases]
		
		self.create_csv(forward_logfc, reverse_logfc)


		return self
		
	def create_csv(self, forward_logfc, reverse_logfc):
		output = []
		for i in range(0,len(forward_logfc)):
			output.append('{} {}\n'.format(forward_logfc[i], reverse_logfc[i]))
		
		with open(self.output_filename, 'w', buffering=1000000) as plotfile:
			plotfile.write(''.join(output))

		return self
		
	def blocks_create(self, logfc_values):
		blocks = []
		inblock = False		
		start = 0
		end = 0
		max_logfc = 0
		
		abs_logfc_values = numpy.absolute(logfc_values)
		for i in range(0,self.genome_length):
			lfc = abs_logfc_values[i]
			if lfc > 0 and not inblock:
				inblock = True
				start = i
				max_logfc = logfc_values[i]
			elif lfc > 0 and inblock:
				if numpy.absolute(max_logfc) < lfc:
					max_logfc = logfc_values[i]
			elif lfc <= 0 and inblock:
				inblock = False
				end = i
				blocks.append(Block(start +1, end, end-start, max_logfc, 'x'))
				max_logfc = 0 
				
		# Check for block at end
		if inblock:
			blocks.append(Block(start +1, len(logfc_values), len(logfc_values)-start, max_logfc, 'x'))
		return blocks	
	
	def genome_wide_logfc(self,logfc_coord_values):
		logfc_to_bases = numpy.zeros(self.genome_length, dtype=float)
		
		# start with the largest signals and overwrite with smaller signals.
		# prevents issues with overlapping blocks/genes
		#sorted_logfc_coord_values = sorted(logfc_coord_values, key = lambda l: (numpy.absolute(l.logfc_value)))
		for l in logfc_coord_values:
			#abs_logfc_value = numpy.absolute(l.logfc_value)
			
			for i in range(l.start -1, l.end):
				logfc_to_bases[i] = l.logfc_value
			
		return self.span_block_gaps(self.filter_out_small_blocks(logfc_to_bases))
			
	def span_block_gaps(self,logfc_to_bases):
		logfc_blocks = self.blocks_create(logfc_to_bases)
		# span blocks if they are close together
		if self.span_gaps > 0:
			for b in logfc_blocks:
				span_index = b.end + (self.window_size * self.span_gaps)
				if span_index >= self.genome_length:
					continue
				
				span_value = numpy.absolute(logfc_to_bases[span_index])
				if span_value >= self.minimum_logfc:
					for a in range(b.end, span_index):
						if numpy.absolute(logfc_to_bases[a]) > self.minimum_logfc:
							continue
						
						if logfc_to_bases[b.end -1 ] < 0:
							logfc_to_bases[a] = -1 * self.minimum_logfc	
						elif logfc_to_bases[b.end -1 ] > 0:
							logfc_to_bases[a] = self.minimum_logfc
					
		return logfc_to_bases
		
	def filter_out_small_blocks(self, logfc_to_bases):
		logfc_blocks = self.blocks_create(logfc_to_bases)
		# filter out small blocks
		for b in logfc_blocks:
			if b.block_length < self.window_size:
				for i in range(b.start -1, b.end):
					logfc_to_bases[i] = 0
					
		return logfc_to_bases
		
	def read_comparison_file(self):
		logfc_coord_values = []
		
		genes_to_features = EMBLReader(self.embl_file).genes_to_features
		
		genes_seen = {}
		with open(self.comparison_filename, newline='') as csvfile:
			comparison_reader = csv.reader(csvfile, delimiter=',')
			for row in comparison_reader:
				logfc = row[3]
				if logfc == 'logFC':
					 continue
					 
				logfc = float(row[3])
				temp_pval = float(row[5])
				temp_logcpm = float(row[4])
				temp_qval = float(row[6])


				
				if not self.report_decreased_insertions and logfc < 0:
					logfc = 0
				
				if numpy.absolute(logfc) < self.minimum_logfc or temp_pval >= self.pvalue:
					logfc = 0
					
				if numpy.absolute(temp_logcpm) < self.minimum_logcpm:
					logfc = 0

				if temp_qval >= self.qvalue:
					logfc = 0


					
				# Get coordinates
				gene_name = row[1]
				
				# else annotation encodes the coordinates
				if gene_name in genes_to_features:
					start = genes_to_features[gene_name].location.start +1
					end = genes_to_features[gene_name].location.end
					logfc_coord_values.append(LogFC(int(start), int(end), logfc))
					# A gene should only be identified once
					genes_seen[gene_name] = 1
				else:
					print("Couldnt find gene coordinates for "+ str(gene_name))
					
			# loop over all the remaining gene and give them a value of zero
			for gene_name,feature in genes_to_features.items():
				if gene_name not in genes_seen:
					start = feature.location.start +1
					end = feature.location.end
					logfc_coord_values.append(LogFC(int(start), int(end), 0))
			
		return logfc_coord_values
		
