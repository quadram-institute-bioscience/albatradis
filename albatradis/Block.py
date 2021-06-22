class Block:
	def __init__(self, start, end, block_length, max_logfc, expression, min_pvalue=1.0, min_qvalue=1.0):
		self.start = start
		self.end = end
		self.max_logfc = max_logfc
		self.min_pvalue = min_pvalue
		self.min_qvalue = min_qvalue
		self.block_length = block_length
		self.expression = expression
		self.direction = 'NA'
		self.genes = []
		self.intergenic = False
		self.num_genes = 0
		self.upstream = 'NA'
		
	def __str__(self):
		return "\t".join(['NA', 'Intergenic', str(self.start), str(self.end), str(self.max_logfc), self.expression, self.direction, self.upstream, str(self.min_pvalue), str(self.min_qvalue)])
