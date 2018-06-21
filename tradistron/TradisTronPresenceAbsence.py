'''Driver class'''
import logging
import os
import sys
import time

from Bio import SeqIO
import csv

class TradisTronPresenceAbsence:
	def __init__(self, options):
		self.logger            = logging.getLogger(__name__)
		self.genereports       = options.genereports
		self.verbose           = options.verbose
		self.emblfile          = options.emblfile
		
		self.genome_length = 0
		
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
			
		if not os.path.exists(self.prefix ):
			os.makedirs(self.prefix )
			
		self.blocks = []
	
	def read_annotation_features(self):
		record =  SeqIO.read(self.emblfile, "embl")
		return [f for f in record.features if f.type not in ['source','gene']]
	
	
	def read_csv(self, filename):
		data = {}
		with open(filename) as csvfile:
			reader = csv.reader(csvfile, delimiter='	')
			data =  [r[0]: r[4] for r in reader]
		return data	
	

	def run(self):
		features = self.read_annotation_features()
		sample_data = {}
		for g in self.genereports:
			sample_data[g] = self.read_csv(g)
		
		gene_names = [f.qualifiers["gene"][0] for f in features]
		print("\t".join(gene_names))
		
		for g in self.genereports:
			row = []
			for gene_name in gene_names:
				if gene_name in g:
					row.append(g[gene_name])
				else:
					row.append(0)
			print("\t".join(row))
		
		return self
