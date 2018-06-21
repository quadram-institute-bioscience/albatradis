import os
import sys
from tradistron.EMBLReader import EMBLReader
from tradistron.GeneReport import GeneReport

class PresenceAbsence:
	def __init__(self, genereports, emblfile, filter_no_data,  verbose ):
		self.genereports       = genereports
		self.verbose           = verbose
		self.emblfile          = emblfile
		self.filter_no_data    = filter_no_data
		
		self.features   = EMBLReader(self.emblfile).read_annotation_features()
		self.gene_names = self.generate_gene_names()
		
	def generate_gene_names(self):
		return [f.qualifiers["gene"][0] for f in self.features if "gene" in f.qualifiers]

	def create_gene_logfc_spreadsheet(self, filename):
		with open(filename, 'w') as fh:
			# Header
			fh.write("\t".join(['Sample'] + self.gene_names) + "\n")

			# Body
			for g in self.genereports:
				genes_to_logfc = GeneReport(g).genes_to_logfc(self.gene_names)
				fh.write("\t".join([g] + genes_to_logfc) + "\n")
				
		return self
