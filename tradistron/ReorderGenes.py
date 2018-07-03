import os
import numpy

from tradistron.GeneToFiles import GeneToFiles

class ReorderGenes:
	def __init__(self, gene_names, genereports, reports_to_gene_logfc):
		self.gene_names            = gene_names
		self.genereports           = genereports
		self.reports_to_gene_logfc = reports_to_gene_logfc
		
		self.genes_to_files = self.create_gene_file_objects()
		
	def create_gene_file_objects(self):
		genes_to_files = {}
		for i,g in enumerate(self.gene_names):
			logfc = []
			for report_file in self.genereports:
				logfc.append(int(self.reports_to_gene_logfc[report_file][i]))
			genes_to_files[g] = GeneToFiles(g, gene_to_files = logfc)	
				
		return genes_to_files
	
	def get_highest_freq(self, genes_to_files):
		sorted_freq = sorted([g.number_of_files for g in genes_to_files.values()], reverse=True)

		return sorted_freq[0]
		
		
			