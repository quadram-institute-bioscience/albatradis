import csv

'''Read in a Gene report CSV file and return the logFC values against a set list of genes'''
class GeneReport:
	def __init__(self, filename):
		self.filename = filename
		self.gene_all_data = self.read_all_data_csv()
		self.gene_data = self.read_csv(self.gene_all_data)
		
	def read_all_data_csv(self):
		all_data = []
		with open(self.filename) as csvfile:
			reader = csv.reader(csvfile, delimiter='	')
			all_data =  [r for r in reader if r[0] != 'Gene']
		return all_data

	def fix_sign_on_logfc(self,gene_all_data):
		for r in gene_all_data:
			r[4] = float(r[4])
			if r[1] == 'upregulated':
				if r[4] < 0.0:
					r[4] *= -1.0
			elif r[1] == 'downregulated':
				if r[4] > 0.0:
					r[4] *= -1.0
			elif r[5] == 'increased_insertions':
				if r[4] < 0.0:
					r[4] *= -1.0
			elif r[5] == 'decreased_insertions':
				if r[4] > 0.0:
					r[4] *= -1.0
		return gene_all_data
		
	def read_csv(self, gene_all_data):
		return {r[0]: r[4] for r in self.fix_sign_on_logfc(gene_all_data)}
		
	def genes_to_logfc(self, gene_names):
		row = []
		for gene_name in gene_names:
			if gene_name in self.gene_data:
				row.append(str(self.gene_data[gene_name]))
			else:
				row.append(str(0))
		return row
		
