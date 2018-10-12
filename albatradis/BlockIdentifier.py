from albatradis.PlotParser import PlotParser
from albatradis.Block import Block
import numpy

class BlockIdentifier:
	
	def __init__(self, combined_mask_file, forward_mask_file, reverse_mask_file, window_size):
		self.combined_mask_file = combined_mask_file
		self.forward_mask_file = forward_mask_file
		self.reverse_mask_file = reverse_mask_file
		self.window_size = window_size
		self.logfc_direction_change = 1
	
	def increased_insertions_blocks(self, masking_plot):
		blocks = []
		inblock = False		
		start = 0
		end = 0
		max_logfc = 0 
		
		for i, mask in enumerate(masking_plot.combined):
			if mask > 0 and not inblock:
				inblock = True
				start = i
				max_logfc = mask
			elif mask > 0 and inblock:
				if max_logfc < mask:
					max_logfc = mask
			elif mask <= 0 and inblock:
				inblock = False
				end = i
				blocks.append(Block(start +1, end, end-start, max_logfc, 'increased_insertions'))
				max_logfc = 0 
				
		# Check for block at end
		if inblock:
			blocks.append(Block(start +1, len(masking_plot.combined), len(masking_plot.combined)-start, max_logfc, 'increased_insertions'))
		return blocks
		
	def decreased_insertions_blocks(self, masking_plot):
		blocks = []
		inblock = False		
		start = 0
		end = 0
		max_logfc = 0
			
		for i, mask in enumerate(masking_plot.combined):
			if mask < 0 and not inblock:
				inblock = True
				start = i
				max_logfc = mask
			elif mask < 0 and inblock:
				if max_logfc > mask:
					max_logfc = mask
			elif mask >= 0 and inblock:
				inblock = False
				end = i
				blocks.append(Block(start +1, end, end-start, max_logfc, 'decreased_insertions' ))
				max_logfc = 0 

		# Check for block at end
		if inblock:
			blocks.append(Block(start +1, len(masking_plot.combined), len(masking_plot.combined)-start, max_logfc, 'decreased_insertions'))
		
		return blocks
		
	def direction_for_block(self, block, forward_masking_plot, reverse_masking_plot):
		forward_max_logfc = 0 
		for i in range(block.start -1, block.end):
			if numpy.absolute(forward_masking_plot.combined[i]) > forward_max_logfc:
				forward_max_logfc = numpy.absolute(forward_masking_plot.combined[i])
				
		reverse_max_logfc = 0 
		for i in range(block.start -1, block.end):
			if numpy.absolute(reverse_masking_plot.combined[i]) > reverse_max_logfc:
				reverse_max_logfc = numpy.absolute(reverse_masking_plot.combined[i])
		
		if forward_max_logfc > reverse_max_logfc:
			if reverse_max_logfc == 0:
				return 'forward'
			elif forward_max_logfc >= reverse_max_logfc + self.logfc_direction_change:
				return 'forward'
			else:
				return 'nodirection'
		else:
			if forward_max_logfc == 0:
				return 'reverse'
			elif reverse_max_logfc >= forward_max_logfc + self.logfc_direction_change:
				return 'reverse'
			else:
				return 'nodirection'

	def block_generator(self):
		masking_plot = PlotParser(self.combined_mask_file)
		blocks = self.increased_insertions_blocks(masking_plot) + self.decreased_insertions_blocks(masking_plot)
		
		'''Filter out blocks which are less than the window size'''
		filtered_blocks = [b for b in blocks if b.block_length >= self.window_size]
		
		forward_masking_plot = PlotParser(self.forward_mask_file)
		reverse_masking_plot = PlotParser(self.reverse_mask_file)
		for b in filtered_blocks:
			b.direction = self.direction_for_block(b, forward_masking_plot, reverse_masking_plot)
		
		return filtered_blocks

# is it beside an essential region
# number of reads in block
# number of insertions
# average reads per insertion
