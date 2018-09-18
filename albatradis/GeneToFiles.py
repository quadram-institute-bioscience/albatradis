import numpy
class GeneToFiles:
	def __init__(self, gene_name, gene_to_files = []):
		self.gene_name = gene_name
		self.gene_to_files = gene_to_files
		self.number_of_files = self.calc_number_of_files()
		
	def calc_number_of_files(self):
		freq = 0
		for g in self.gene_to_files:
			if numpy.absolute(g) > 0:
				freq += 1
		return freq
		
	def __str__(self):
		return "\t".join([self.gene_name, str(self.number_of_files)])