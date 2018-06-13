class Block:
	def __init__(self, start, end, block_length, max_logfc, expression, direction ):
		self.start = start
		self.end = end
		self.max_logfc = max_logfc
		self.block_length = block_length
		self.expression = expression
		self.direction = direction	