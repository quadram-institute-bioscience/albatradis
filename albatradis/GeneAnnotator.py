from Bio import SeqIO
import re
from albatradis.Gene import Gene
from albatradis.EMBLReader import EMBLReader

class GeneAnnotator:
	def __init__(self, annotation_file, blocks):
		self.annotation_file = annotation_file
		self.blocks = self.sort_blocks_by_start_coord(blocks)
		self.knockout_proportion_start = 0.5
		self.increased_expression_proportion_end= 0.3
		
		self.embl_reader =  EMBLReader(self.annotation_file)
		self.features = self.embl_reader.features
		
	def sort_blocks_by_start_coord(self, blocks):
		return sorted((b for b in blocks ), key=lambda x: x.start)
		
	def annotate_genes(self):
		genes = []
		for gene_number, f in enumerate(self.features):
			overlapping_blocks = self.blocks_overlapping_feature(f) 

			if len(overlapping_blocks) == 0:
				# no hits to any blocks so move to next feature
				continue

			g = Gene(f, overlapping_blocks)
			 
			# only consider block at a time
			for b in overlapping_blocks:
				g.upstream.append(self.find_upstream_gene(b,gene_number))
				if self.is_feature_contained_within_block(b, f):
					g.categories.append('knockout')
				elif self.is_block_near_end_of_feature(b, f):
					if b.max_logfc > 0 :
						g.categories.append('increased_mutants_at_end_of_gene')
					else:
						g.categories.append('decreased_mutants_at_end_of_gene')
				elif self.is_block_near_start_of_feature(b,f):
					if b.max_logfc > 0 :
						g.categories.append('increased_mutants_at_start_of_gene')
					else:
						g.categories.append('decreased_mutants_at_start_of_gene')
			
			if len(g.categories) == 0:
				p = self.proportion_blocks_overlap_with_gene(f, overlapping_blocks)
				if p > 0.9:
					g.categories.append('over_90_perc_inactivation')
				elif p > 0.8:
					g.categories.append('over_80_perc_inactivation')	
				elif p > 0.7:
					g.categories.append('over_70_perc_inactivation')
				elif p > 0.6:
					g.categories.append('over_60_perc_inactivation')
				elif p > 0.5:
					g.categories.append('over_50_perc_inactivation')	
			
			if len(g.categories) == 0:
				g.categories.append('unclassified')
			genes.append(g)

		# intergenic test
		intergenic_blocks = [block for block in self.blocks if block.num_genes == 0]
		for block in intergenic_blocks:
			block.upstream  = self.find_nearest_upstream_gene(block)
			block.intergenic = True

		reannotate_with_5_3_prime = self.reannotate_5_3_prime(genes)

		return reannotate_with_5_3_prime 
		
	def feature_to_gene_name(self, feature):
		gene_name_val = str(feature.location.start) + "_" + str(feature.location.end)
		if "gene" in feature.qualifiers:
			gene_name_val = feature.qualifiers["gene"][0]
		return gene_name_val
		
	def find_nearest_upstream_gene(self, block):
		if block.direction == 'forward':
			for f in self.features:
				if f.location.start -1 > block.end and f.strand == 1:
					return self.feature_to_gene_name(f)
		elif block.direction == 'reverse':
			for f in reversed(self.features):
				if f.location.end < block.start and f.strand == -1:
					return self.feature_to_gene_name(f)
		return "NA"
		
	def find_upstream_gene(self, block, gene_number):
		if block.direction  == 'reverse':
			for upstream_i in range(gene_number -1, 0, -1):
				if self.features[upstream_i].strand == -1 and block.start -1  > self.features[upstream_i].location.end:
					return self.feature_to_gene_name(self.features[upstream_i])
		elif block.direction  == 'forward':
			for upstream_i in range(gene_number +1, len(self.features)):
				if self.features[upstream_i].strand == 1 and block.end < self.features[upstream_i].location.start:
					return self.feature_to_gene_name(self.features[upstream_i])
		return "NA"
		
	def proportion_blocks_overlap_with_gene(self,gene, blocks):
		base_coverage = 0
		for b in blocks:
			for b_index in range (b.start -1, b.end):
				if b_index >= gene.location.start and b_index < gene.location.end:
					base_coverage += 1
				
		gene_length = gene.location.end - gene.location.start
		return base_coverage/gene_length
			
		
	def blocks_overlapping_feature(self, feature):
		overlapping_blocks = []
		
		for block in self.blocks:

			if (block.start -1) > feature.location.end  or feature.location.start > block.end:
				continue
				
			# genes are big so you are bound to hit one. Smallest in ecoli is 45bp so half it.
			for i in range(block.start -1 , block.end, 22):
				if i in feature:
					overlapping_blocks.append(block)
					block.num_genes += 1
					break
		return overlapping_blocks
			
	def is_feature_contained_within_block(self, block, feature):
		if feature.location.start >= block.start -1 and feature.location.end <= block.end:
			return True
		return False
		
	def is_block_contained_within_feature(self, block, feature):
		if block.start -1 >= feature.location.start and block.end <= feature.location.end:
			return True
		return False
		
	def is_block_overlapping_feature_on_right(self, block, feature):
		if block.start -1  < feature.location.end and block.start -1 > feature.location.start and block.end > feature.location.end:
			return True
		return False
		
	def is_block_overlapping_feature_on_left(self, block, feature):	
		if block.end < feature.location.end and block.end > feature.location.start and block.start -1 < feature.location.start:
			return True
		return False
		
	def is_block_near_end_of_feature(self, block, feature):
		# forward 
		if feature.strand == 1 and block.direction in ['reverse', 'nodirection']:
			knock_out_start = feature.location.start + int(self.knockout_proportion_start * len(feature))
			if block.start -1  >= knock_out_start and block.start -1  < feature.location.end:
				return True
		# reverse
		if feature.strand == -1 and block.direction in ['forward', 'nodirection']:		
			knock_out_end = feature.location.end - int(self.knockout_proportion_start * len(feature))
			if block.end <= knock_out_end and block.end  > feature.location.start:
				return True
		
		return False
			
	def is_block_near_start_of_feature(self, block, feature):
		# forward
		if feature.strand == 1 and block.direction in ['forward','nodirection']:
			expression_end = feature.location.start +  int(self.increased_expression_proportion_end * len(feature))
			if block.end <= expression_end and block.end > feature.location.start:
				return True
		# reverse
		if feature.strand == -1 and block.direction in ['reverse', 'nodirection']:
			expression_end = feature.location.start +  int(self.increased_expression_proportion_end * len(feature))
			if block.end <= expression_end and block.end > feature.location.start:
				return True
				
		return False
		
		
	def regulation(self,strand, prime, directions):
		if prime == '5' and strand == 1 and 'forward' in directions:
			return 'upregulated'
		elif prime == '5' and strand == -1 and 'reverse' in directions:
			return 'upregulated'
		elif prime == '3' and strand == 1 and 'reverse' in directions:
			return 'downregulated'
		elif prime == '3' and strand == -1 and 'forward' in directions:
			return 'downregulated'
		else:
			return None
		
	def reannotate_5_3_prime(self,genes):
		name_to_genes = { g.gene_name: g for g in genes}
		
		filtered_names_to_genes = {}

		for name, gene in name_to_genes.items():
			directions = list(set([b.direction for b in gene.blocks]))
			
			if 'nodirection' in directions:
				continue	
			# 5 prime
			res = re.search("^(.+)__([35])prime$", name)
			if res:
				found_gene_name = res.group(1)
				prime_end = res.group(2)
				if found_gene_name not in name_to_genes:
					filtered_names_to_genes[found_gene_name] = Gene(self.embl_reader.genes_to_features[found_gene_name], [])
					filtered_names_to_genes[found_gene_name].blocks = gene.blocks
					
					
					regulation_category = self.regulation(filtered_names_to_genes[found_gene_name].feature.strand,prime_end, directions)
					if regulation_category:
						filtered_names_to_genes[found_gene_name].categories.append(regulation_category)
					else:
						del filtered_names_to_genes[found_gene_name]
					
				else:
					filtered_names_to_genes[found_gene_name] = name_to_genes[found_gene_name]
					regulation_category = self.regulation(filtered_names_to_genes[found_gene_name].feature.strand,prime_end, directions)
					if regulation_category:
						filtered_names_to_genes[found_gene_name].categories.append(regulation_category)

		# carry over non prime genes
		for name, gene in name_to_genes.items():
			res = re.search("^(.+)__[35]prime$", name)
			if not res:
				if name not in filtered_names_to_genes:
					filtered_names_to_genes[name] = gene
					
		filtered_names_to_genes = self.fix_pos_neg(filtered_names_to_genes)

		return [g for g in filtered_names_to_genes.values()]
		
	def fix_pos_neg(self, filtered_names_to_genes):
		for name, gene in filtered_names_to_genes.items():
			# If there are increased insertions make it pos, decreased make it neg
			increased = [True for g in filtered_names_to_genes[name].blocks if g.expression == 'increased_insertions']
			for b in filtered_names_to_genes[name].blocks:
				if b.max_logfc < 0 and True in increased:
					b.max_logfc *= -1
					
			decreased = [True for g in filtered_names_to_genes[name].blocks if g.expression == 'decreased_insertions']
			for b in filtered_names_to_genes[name].blocks:
				if b.max_logfc > 0 and True in decreased:
					b.max_logfc *= -1
			
			# If there is downregulation make it negative
			downreg = [True for g in filtered_names_to_genes[name].categories if 'downregulated' == g]
			for b in filtered_names_to_genes[name].blocks:
				if b.max_logfc >= 0 and True in downreg:
					b.max_logfc *= -1
					
			# If there is upregulation make it positive
			upreg = [True for g in filtered_names_to_genes[name].categories if 'upregulated' == g]
			for b in filtered_names_to_genes[name].blocks:
				if b.max_logfc < 0 and True in upreg:
					b.max_logfc *= -1
		return filtered_names_to_genes
