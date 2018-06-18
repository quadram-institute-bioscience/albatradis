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
			output_strs.append("\t".join([str(gene.id), str(gene.category),str(self.start), str(self.end), str(self.block_length), str(self.max_logfc),  str(self.expression), str(self.direction), str(self.category)] ))
		return join("\n", output_strs)
		
	def header(self):
		"\t".join(['Start', 'End','Length', 'MaxLogFC', 'Expression', 'Direction', 'Category'])
		