class Block:
	def __init__(self, start, end, block_length, max_logfc, expression ):
		self.start = start
		self.end = end
		self.max_logfc = max_logfc
		self.block_length = block_length
		self.expression = expression
		self.direction = 'unknown'
		self.genes = []
		
	def __str__(self):
		output_strs = []
		for gene in self.genes:
			
			gene_name = 'unknown'
			if "gene" in gene.feature.qualifiers:
				gene_name = gene.feature.qualifiers["gene"][0]
			
			output_strs.append("\t".join([str(gene_name), str(gene.category),str(self.start), str(self.end), str(self.block_length), str(self.max_logfc),  str(self.expression), str(self.direction)] ))
		else:
			output_strs.append("\t".join(['unknown', '',str(self.start), str(self.end), str(self.block_length), str(self.max_logfc),  str(self.expression), str(self.direction)] ))
		return "\n".join( output_strs)
		
	def header(self):
		return "\t".join(['Gene', 'Category', 'Start', 'End','Length', 'MaxLogFC', 'Expression', 'Direction', 'Category'])
		