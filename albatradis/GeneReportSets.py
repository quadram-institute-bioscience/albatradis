from albatradis.GeneReport import GeneReport
from albatradis.Gene import Gene
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
	
	# Find the genes that are present in each file
	def intersection(self):
		combined = {}
		candidate_genes = {}
		first_filename = self.filenames[0]
		# get all the gene names that are present in the 1st file
		for row in self.gene_reports[first_filename].gene_all_data:
			gene = self.row_to_gene_name(row)
			combined[gene] = row
			candidate_genes[gene] = 1
		
		# Count the number of times each of the 1st set of genes occurs
		for f in self.filenames:
			if f == first_filename:
				continue
			for row in self.gene_reports[f].gene_all_data:
				gene = self.row_to_gene_name(row)
					
				if gene in combined:
					candidate_genes[gene] += 1
		
		# Extract only the genes which occur in every file and map to their rows in the spreadsheet
		genes_found_in_all  = [genename for genename,count in candidate_genes.items() if count ==  len(self.filenames)]
		filtered_gene_rows = {}
		for g in genes_found_in_all:
			filtered_gene_rows[g] = combined[g]
		return filtered_gene_rows
	
	def row_to_gene_name(self,row):
		gene = row[0]
		# use the start and end coords for unnamed features
		if gene == 'unknown' or gene == 'NA':
			gene = str(row[2]) + "_" + str(row[3])
		return gene
					
	# When merging, use the first row for a gene. This can have unintended consequences (like an increase in insertions in one exp, and a decrease in insertions in another)
	def union(self):
		combined = {}
		for f in self.filenames:
			for row in self.gene_reports[f].gene_all_data:
				gene = self.row_to_gene_name(row)
					
				if gene not in combined:
					combined[gene] = row
					
		return combined
					
	def write_union_file(self):
		union_filename = os.path.join(self.prefix, "union_gene_report.csv")
		
		with open(union_filename, 'w') as bf:
			bf.write(str(Gene(None,[]).header())+"\n")
			for row in sorted(self.union().values(), key=lambda x: int(x[2])):
				bf.write( "\t".join([str(i) for i in row]) +"\n")
				
		return self

	def write_intersection_file(self):
		intersection_filename = os.path.join(self.prefix, "intersection_gene_report.csv")
		
		if len(self.intersection()) > 0:
			with open(intersection_filename, 'w') as bf:
				bf.write(str(Gene(None,[]).header())+"\n")
				for row in sorted(self.intersection().values(), key=lambda x: int(x[2])):
					bf.write( "\t".join([str(i) for i in row]) +"\n")
		else:
			print("No intersecting genes")
			

				
		return self
