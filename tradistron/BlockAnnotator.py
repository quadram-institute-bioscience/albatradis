from Bio import SeqIO
from tradistron.Gene import Gene

class BlockAnnotator:
	def __init__(self, annotation_file, blocks):
		self.annotation_file = annotation_file
		self.blocks = blocks
		self.knockout_proportion_start = 0.5
		self.increased_expression_proportion_end= 0.3
		
		self.features = self.read_annotation_features()
	
	def annotate_blocks(self):
		for b in self.blocks:
			overlapping_feature = self.features_overlapping_block(b) 
			
			# only consider 1 gene at a time
			for f in overlapping_feature:
				if self.is_feature_contained_within_block(b, f):
					b.genes.append(Gene(f, 'total_inactivation'))
					continue
				elif self.is_block_near_end_of_feature(b, f):
					if b.max_logfc > 0 :
						b.genes.append(Gene(f, 'increased_mutants_at_end_of_gene'))
					else:
						b.genes.append(Gene(f, 'decreased_mutants_at_end_of_gene'))
					continue
				elif self.is_block_near_start_of_feature(b,f):
					if b.max_logfc > 0 :
						b.genes.append(Gene(f, 'increased_mutants_at_start_of_gene'))
					else:
						b.genes.append(Gene(f, 'decreased_mutants_at_start_of_gene'))
					continue
					
				
		return self 
		
	'''Assumption is that there is only 1 sequence in the annotation file'''
	def read_annotation_features(self):
		record =  SeqIO.read(self.annotation_file, "embl")
		
		return record.features
		
	def features_overlapping_block(self,block):
		overlapping_features = []
		for f in self.features:
			if f.type in ['source','gene']:
				continue
			
			for i in range(block.start -1 , block.end):
				if i in f:
					overlapping_features.append(f)
					break
		return overlapping_features
			
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
		
		
		
		
		
		
		