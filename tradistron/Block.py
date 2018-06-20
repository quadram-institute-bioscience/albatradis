class Block:
	def __init__(self, start, end, block_length, max_logfc, expression ):
		self.start = start
		self.end = end
		self.max_logfc = max_logfc
		self.block_length = block_length
		self.expression = expression
		self.direction = 'NA'
		self.genes = []
		self.intergenic = False
		self.num_genes = 0
		
	def __str__(self):
		return "\t".join(['NA', 'Intergenic', str(self.start), str(self.end), str(self.max_logfc), self.expression, self.direction])
