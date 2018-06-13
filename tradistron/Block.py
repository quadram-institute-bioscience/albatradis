class Block:
	def __init__(self, start, end, block_length, max_logfc, expression ):
		self.start = start
		self.end = end
		self.max_logfc = max_logfc
		self.block_length = block_length
		self.expression = expression
		self.direction = 'unknown'
		
	def __str__(self):
		return "\t".join([str(self.start), str(self.end), str(self.block_length), str(self.max_logfc),  str(self.expression), str(self.direction)] )
		
	def header(self):
		"\t".join(['Start', 'End','Length', 'MaxLogFC', 'Expression', 'Direction'])
		