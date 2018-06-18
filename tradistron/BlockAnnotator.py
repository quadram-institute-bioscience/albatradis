from Bio import SeqIO
from tradistron.Gene import Gene

class BlockAnnotator:
	def __init__(self, annotation_file, blocks):
		self.annotation_file = annotation_file
		self.blocks = blocks
		
		self.features = self.read_annotation_features()
	
	def annotate_blocks(self):
		for b in self.blocks:
			overlapping_feature = self.features_overlapping_block(b) 
			for f in overlapping_feature:
				if self.is_feature_contained_within_block(b, f):
					b.genes.append(Gene(f, 'total_inactivation'))
					continue
		
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
		
		