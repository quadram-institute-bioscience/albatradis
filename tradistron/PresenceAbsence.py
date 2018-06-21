import os
import sys
import numpy
from tradistron.EMBLReader import EMBLReader
from tradistron.GeneReport import GeneReport

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt


class PresenceAbsence:
	def __init__(self, genereports, emblfile, verbose, prefix):
		self.genereports       = genereports
		self.verbose           = verbose
		self.emblfile          = emblfile
		self.prefix            = prefix
		
		self.run()
		
	def run(self):
		self.features   = EMBLReader(self.emblfile).read_annotation_features()
		self.gene_names = self.generate_gene_names()
		self.reports_to_gene_logfc = self.create_reports_to_gene_logfc(self.gene_names)
		
		self.filtered_gene_names = self.filter_genes_with_no_changes()
		self.filtered_reports_to_gene_logfc = self.create_reports_to_gene_logfc(self.filtered_gene_names)
		
	def generate_gene_names(self):
		return [f.qualifiers["gene"][0] for f in self.features if "gene" in f.qualifiers]

	def create_reports_to_gene_logfc(self, gene_names):
		reports_to_gl = {}
		for g in self.genereports:
			reports_to_gl[g] = GeneReport(g).genes_to_logfc(gene_names)
		return 	reports_to_gl

	def create_output_files(self):
		self.create_gene_logfc_spreadsheet(os.path.join(self.prefix, 'all_logfc.csv'), self.gene_names, self.reports_to_gene_logfc)
		self.create_gene_logfc_spreadsheet(os.path.join(self.prefix, 'filtered_logfc.csv'), self.filtered_gene_names, self.filtered_reports_to_gene_logfc)
		
		self.plot_distance_matrix(os.path.join(self.prefix, 'distance_matrix_dendrogram.png'))
		return self

	def create_gene_logfc_spreadsheet(self, filename, gene_names, reports_to_gene_logfc):
		with open(filename, 'w') as fh:
			# Header
			fh.write("\t".join(['Sample'] + gene_names) + "\n")

			# Body
			for g in self.genereports:
				fh.write("\t".join([g] + reports_to_gene_logfc[g]) + "\n")
				
		return self
		
	def filter_genes_with_no_changes(self):
		genes_with_changes = []
		for i,g in enumerate(self.gene_names):
			for report_file in self.genereports:
				if numpy.absolute(int(self.reports_to_gene_logfc[report_file][i])) > 0:
					genes_with_changes.append(g)
					break
		return genes_with_changes

	# assumption is that any changes (+-) on a gene is counted the same.
	def pair_wise_distance(self, file_a, file_b):
		distance = 0 
		for i in range(len(self.filtered_gene_names)):
			a_abs = numpy.absolute(int(self.filtered_reports_to_gene_logfc[file_a][i]))
			b_abs = numpy.absolute(int(self.filtered_reports_to_gene_logfc[file_b][i]))
			
			if (a_abs == 0 and b_abs == 0) or (a_abs > 0 and b_abs > 0):
				continue
			else:
				distance += 1 
		return distance

	def distance_matrix(self):
		distances = []
		for a in self.genereports:
			dist_row = []
			for b in self.genereports:
				dist_row.append(self.pair_wise_distance(a,b))
			distances.append(dist_row)
		return distances
		
	def plot_distance_matrix(self, outputfile):
		mat = numpy.array(self.distance_matrix())
		dists = squareform(mat)
		linkage_matrix = linkage(dists, "single")
		dendrogram(linkage_matrix, labels=self.genereports )
		plt.title("Distance matrix")
		plt.savefig(outputfile)
		
		return self
			
				