import numpy
class Gene:
	def __init__(self, feature, blocks):
		self.feature  = feature
		self.blocks = blocks
		self.categories = []
		self.upstream = []
		self.five_prime = None
		self.three_prime = None
		
		self.gene_name = self.calc_gene_name()

	def calc_gene_name(self):
		if self.feature:
			gene_name_val = str(self.feature.location.start) + "_" + str(self.feature.location.end)
			if "gene" in self.feature.qualifiers:
				gene_name_val = self.feature.qualifiers["gene"][0]
				return gene_name_val
		else:
			return "unknown"
		
	def category(self):
		if 'knockout' in self.categories:
			return 'knockout'
		elif 'upregulated' in self.categories:
			return 'upregulated'
		elif 'downregulated' in self.categories:
			return 'downregulated'
		else:
			return "/".join(list(set(self.categories)))
			
	def upstream_gene(self):
		return "/".join(list(set(self.upstream)))

	def max_logfc_from_blocks(self):
		if self.blocks:
			all_logfc = [b.max_logfc for b in self.blocks]
			highest_logfc = numpy.max(numpy.absolute(all_logfc))
			for a in all_logfc:
				if a == highest_logfc:
					return highest_logfc
				elif a == highest_logfc*-1:
					return highest_logfc*-1
		else:
			return 0

	def expression_from_blocks(self):
		if self.blocks:
			return "/".join(list(set([b.expression for b in self.blocks])))
		else:
			return ""
		
	def direction_from_blocks(self):
		if self.blocks:
			return "/".join(list(set([b.direction for b in self.blocks])))
		else:
			return ""

	def __str__(self):
		return "\t".join([str(self.gene_name), str(self.category()), str(self.feature.location.start), str(self.feature.location.end), str(self.max_logfc_from_blocks()),  str(self.expression_from_blocks()), str(self.direction_from_blocks()), str(self.upstream_gene())] )
		
	def header(self):
		return "\t".join(['Gene', 'Category', 'Start', 'End', 'MaxLogFC', 'Expression', 'Direction', 'Upstream'])
		
