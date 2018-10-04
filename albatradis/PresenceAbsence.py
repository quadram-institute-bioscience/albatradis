import os
import sys
import numpy
import dendropy

from dendropy.utility.textprocessing import StringIO
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from graphviz import Digraph

from albatradis.EMBLReader import EMBLReader
from albatradis.GeneReport import GeneReport
from albatradis.HeatMap import HeatMap
from albatradis.GeneToFiles import GeneToFiles
from albatradis.ReorderGenes import ReorderGenes

class PresenceAbsence:
	def __init__(self, genereports, emblfile, verbose, prefix):
		self.genereports       = genereports
		self.verbose           = verbose
		self.emblfile          = emblfile
		self.prefix            = prefix
		
		self.run()
		
	def run(self):
		self.features   = EMBLReader(self.emblfile).features
		self.gene_names = self.generate_gene_names()
		self.reports_to_gene_logfc = self.create_reports_to_gene_logfc(self.gene_names)
		
		self.filtered_gene_names = ReorderGenes(self.gene_names, self.genereports, self.reports_to_gene_logfc).reorder_genes()
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
		
		self.create_dot_graph_genes(os.path.join(self.prefix, 'logfc.dot'), self.filtered_gene_names, self.filtered_reports_to_gene_logfc)
		
		self.plot_distance_matrix(os.path.join(self.prefix, 'distance_matrix_dendrogram.png'))
		self.create_nj_newick(os.path.join(self.prefix, 'nj_newick_tree.tre'))
		
		HeatMap(self.reports_to_gene_logfc, self.gene_names,os.path.join(self.prefix, 'full_heatmap.png')).create_heat_map()
		HeatMap(self.filtered_reports_to_gene_logfc, self.filtered_gene_names, os.path.join(self.prefix, 'filtered_heatmap.png')).create_heat_map()
		
		return self
		
		
	def create_dot_graph_genes(self, filename, gene_names, reports_to_gene_logfc):
		dot = Digraph(comment='Gene Graph')
		
		# Antibiotics
		for g in self.genereports:
			dot.node(g, g)
		
		# gene names
		for g in gene_names:
			dot.node(g, g)
			
		for anti in self.genereports:
			for i,logfc in enumerate(reports_to_gene_logfc[anti]):
				if numpy.absolute(int(logfc)) > 0:
					dot.edge(anti, gene_names[i])
				
		with open(filename, 'w') as fh:
			fh.write(dot.source)
			
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
		
		gene_to_freq = {}
		for i,g in enumerate(self.gene_names):
			for report_file in self.genereports:
				cell_logfc = numpy.absolute(int(self.reports_to_gene_logfc[report_file][i]))
				if cell_logfc > 0:
					if g in gene_to_freq:
						gene_to_freq[g] += 1 
					else:
						gene_to_freq[g] = 1
		
		# sort by value descending, then by gene name asc	
		sorted_genes = sorted(gene_to_freq.items(), key=lambda x: (-x[1], x[0]))
		gene_name_index = { g:i for i,g in enumerate(self.gene_names) }
		sorted_gene_index = []
		for s in sorted_genes:
			ordered_gene = s[0]
			sorted_gene_index.append(gene_name_index[ordered_gene])
			
		for i in sorted_gene_index:
			g = self.gene_names[i]
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
		
		
	def nj_distance_matrix_str(self):
		output = "\t".join(['.'] + self.genereports) + "\n"
		
		dm = self.distance_matrix()
		for i,g in enumerate(self.genereports):
			
			output += "\t".join([g] + [str(d) for d in dm[i]]) + "\n"
			
		return output
			
		
	def create_nj_newick(self,outputfile):
		pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
		       src=StringIO(self.nj_distance_matrix_str()),
		        delimiter="\t")
		nj_tree = pdm.nj_tree()
		
		with open(outputfile, 'w') as fh:
			fh.write(nj_tree.as_string("newick"))
			
		return self
			
				
