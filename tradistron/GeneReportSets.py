from tradistron.GeneReport import GeneReport
from tradistron.Gene import Gene
import os


'''Take in 2 or more gene report spreadsheets and output the union and intersection, ...'''
class GeneReportSets:
	def __init__(self, filenames, prefix):
		self.filenames = filenames
		self.prefix = prefix
		
		if not os.path.exists(self.prefix ):
			os.makedirs(self.prefix )
			
		self.gene_reports = self.parse_gene_reports()
		
	def parse_gene_reports(self):
		return { filename: GeneReport(filename) for filename in self.filenames }
		
	# When merging, use the first row for a gene. This can have unintended consequences (like an increase in insertions in one exp, and a decrease in insertions in another)
	def union(self):
		combined = {}
		for f in self.filenames:
			for row in self.gene_reports[f].gene_all_data:
				gene = row[0]
				# use the start and end coords for unnamed features
				if gene == 'unknown' or gene == 'NA':
					gene = str(row[2]) + "_" + str(row[3])
					
				if gene not in combined:
					# normal gene name
					combined[gene] = row
					
		return combined
					
	def write_union_file(self):
		union_filename = os.path.join(self.prefix, "union_gene_report.csv")
		
		with open(union_filename, 'w') as bf:
			bf.write(str(Gene(None,[]).header())+"\n")
			for row in sorted(self.union().values(), key=lambda x: int(x[2])):
				bf.write( "\t".join([str(i) for i in row]) +"\n")
				
		return self
