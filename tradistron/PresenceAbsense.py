import os
import sys
from tradistron.EMBLReader import EMBLReader
from tradistron.GeneReport import GeneReport

class PresenceAbsence:
	def __init__(self, genereports, emblfile, filter_no_data,verbose ):
		self.genereports       = genereports
		self.verbose           = verbose
		self.emblfile          = emblfile
		self.filter_no_data = filter_no_data

	def create_gene_logfc_spreadsheet(self):
		features = EMBLReader(self.emblfile).read_annotation_features()
		gene_names = [f.qualifiers["gene"][0] for f in features if "gene" in f.qualifiers]
		print("\t".join(['Sample'] + gene_names))
		
		for g in self.genereports:
			genes_to_logfc = GeneReport(g).genes_to_logfc(gene_names)
			print("\t".join([g] + genes_to_logfc))
		
		return self
