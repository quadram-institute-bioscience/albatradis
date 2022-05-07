from albatradis.GeneToFiles import GeneToFiles


class ReorderGenes:
	def __init__(self, gene_names, gene_reports, gene_list):
		self.gene_names = gene_names
		self.gene_reports = gene_reports
		self.gene_list = gene_list
		
		self.genes_to_files = self.create_gene_file_objects()
		
	def create_gene_file_objects(self):
		genes_to_files = {}
		for i, g in enumerate(self.gene_names):
			logfc = []
			for report_file in self.gene_reports:
				logfc.append(self.gene_list[report_file][i].max_logfc if self.gene_list[report_file][i] is not None else 0.0)

			gf = GeneToFiles(g, gene_to_files=logfc)
			if gf.number_of_files > 0:
				genes_to_files[g] = gf
				
		return genes_to_files
	
	def get_highest_freq(self, genes_to_files):
		if len(genes_to_files) == 0:
			return None
		
		sorted_freq = sorted(genes_to_files.items(), key=lambda kv: kv[1].number_of_files, reverse=True)
		return sorted_freq[0][0]
		
	def files_in_common(self, first_gene, second_gene):
		common = 0 
		for i in range(len(first_gene.gene_to_files)):
			if first_gene.gene_to_files[i] != 0 and second_gene.gene_to_files[i] != 0:
				common += 1
		return common
		
	def find_closest_gene(self, gene_name, genes_to_files):
		genes_to_num_common = {}
		for gf in genes_to_files.values():
			if gene_name == gf.gene_name:
				continue
			genes_to_num_common[gf.gene_name] = self.files_in_common(genes_to_files[gene_name], gf)
		
		sorted_genes_by_common = sorted(genes_to_num_common.items(), key=lambda kv: kv[1], reverse=True)	
		num_in_common = sorted_genes_by_common[0][1]
		
		all_genes_with_highest_num_in_common = { k:v for k,v in genes_to_num_common.items() if v == num_in_common }
		gene_names_to_num_files = { g: genes_to_files[g].number_of_files for g in all_genes_with_highest_num_in_common.keys() }
		
		sorted_gene_names_to_num_files = sorted(genes_to_num_common.items(), key=lambda kv: kv[1], reverse=True)
		return sorted_gene_names_to_num_files[0][0]
		
		
	def reorder_genes(self):
		if len(self.genes_to_files) == 0:
			return []
		highest_freq_gene_name = self.get_highest_freq(self.genes_to_files)
		reordered_genes = [] 

		while len(self.genes_to_files) > 0:	
			if len(self.genes_to_files) == 1:
				reordered_genes.append(self.genes_to_files.pop(highest_freq_gene_name))
				break
			
			closest_gene = self.find_closest_gene(highest_freq_gene_name, self.genes_to_files)
			
			reordered_genes.append(self.genes_to_files.pop(highest_freq_gene_name))
			highest_freq_gene_name = closest_gene
			if len(self.genes_to_files) == 0:
				return []
			
		return [r.gene_name for r in reordered_genes]

