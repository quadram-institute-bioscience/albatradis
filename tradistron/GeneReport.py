import csv

'''Read in a Gene report CSV file and return the logFC values against a set list of genes'''
class GeneReport:
	def __init__(self, filename):
		self.filename = filename
		self.gene_data = self.read_csv()
		
	def read_csv(self):
		data = {}
		with open(self.filename) as csvfile:
			reader = csv.reader(csvfile, delimiter='	')
			data =  {r[0]: r[4] for r in reader}
		return data
		
	def genes_to_logfc(self, gene_names):
		row = []
		for gene_name in gene_names:
			if gene_name in self.gene_data:
				row.append(str(self.gene_data[gene_name]))
			else:
				row.append(str(0))
		return row
		