import csv
from albatradis.Gene import Gene

'''Read in a Gene report CSV file and return the logFC values against a set list of genes'''
class GeneReport:
	def __init__(self, filename):
		self.filename = filename
		self.gene_all_data = self.read_all_data_csv()
		self.gene_data = self.read_csv(self.gene_all_data)
		
	def read_all_data_csv(self):
		with open(self.filename) as genereportfile:
			all_data = [Gene.parseLine(r.strip()) for r in genereportfile if not r.strip().startswith('Gene')]
		return all_data

	def fix_sign_on_logfc(self, gene_all_data):
		for r in gene_all_data:
			if r.categories[0] == 'upregulated':
				if r.max_logfc < 0.0:
					r.max_logfc *= -1.0
			elif r.categories[0] == 'downregulated':
				if r.max_logfc > 0.0:
					r.max_logfc *= -1.0
			elif r.expression == 'increased_insertions':
				if r.max_logfc < 0.0:
					r.max_logfc *= -1.0
			elif r.expression == 'decreased_insertions':
				if r.max_logfc > 0.0:
					r.max_logfc *= -1.0
		return gene_all_data
		
	def read_csv(self, gene_all_data):
		return {r.gene_name: r for r in self.fix_sign_on_logfc(gene_all_data)}

	def filtered_genes(self, gene_names):
		row = []
		for gene_name in gene_names:
			if gene_name in self.gene_data:
				row.append(self.gene_data[gene_name])
			else:
				row.append(None)
		return row

	def genes_to_logfc(self, gene_names):
		row = []
		for gene_name in gene_names:
			if gene_name in self.gene_data:
				row.append(str(self.gene_data[gene_name].max_logfc))
			else:
				row.append(str(0.0))
		return row

	def genes_to_qvals(self, gene_names):
		row = []
		for gene_name in gene_names:
			if gene_name in self.gene_all_data:
				row.append(str(self.gene_all_data[gene_name].min_qvalue))
			else:
				row.append(str(1.0))
		return row
