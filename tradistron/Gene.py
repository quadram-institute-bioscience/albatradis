import numpy
class Gene:
	def __init__(self, feature, blocks):
		self.feature  = feature
		self.blocks = blocks
		self.categories = []

	def gene_name(self):
		gene_name_val = 'unknown'
		if "gene" in self.feature.qualifiers:
			gene_name_val = self.feature.qualifiers["gene"][0]
		return gene_name_val
		
	def category(self):
		if 'total_inactivation' in self.categories:
			return 'total_inactivation'
		else:
			return "/".join(self.categories)

	def max_logfc_from_blocks(self):
		highest_logfc = 0
		for b in self.blocks:
			if numpy.absolute(b.max_logfc) > numpy.absolute(highest_logfc):
				highest_logfc = b.max_logfc
		return highest_logfc
		
	def expression_from_blocks(self):
		return "/".join(list(set([b.expression for b in self.blocks])))
		
	def direction_from_blocks(self):
		return "/".join(list(set([b.direction for b in self.blocks])))

	def __str__(self):
		return "\t".join([str(self.gene_name()), str(self.category()), str(self.feature.location.start), str(self.feature.location.end), str(self.max_logfc_from_blocks()),  str(self.expression_from_blocks()), str(self.direction_from_blocks())] )
		
	def header(self):
		return "\t".join(['Gene', 'Category', 'Start', 'End', 'MaxLogFC', 'Expression', 'Direction'])
		