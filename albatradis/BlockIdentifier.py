import numpy

from albatradis.Block import Block
from albatradis.PlotParser import PlotParser, ScoreParser


class BlockIdentifier:

    def __init__(self, combined_mask_file, forward_mask_file, reverse_mask_file, combined_score_file, window_size):
        self.combined_mask_file = combined_mask_file
        self.forward_mask_file = forward_mask_file
        self.reverse_mask_file = reverse_mask_file
        self.combined_score_file = combined_score_file
        self.window_size = window_size
        self.logfc_direction_change = 1

    def increased_insertions_blocks(self, masking_plot_combined, score_plot_combined):
        blocks = []
        inblock = False
        start = 0
        end = 0
        max_logfc = 0.0
        min_pval = 1.0
        min_qval = 1.0

        for i, mask in enumerate(masking_plot_combined):
            if mask > 0.0 and not inblock:
                inblock = True
                start = i
                max_logfc = mask
                min_pval = score_plot_combined.pvals[i]
                min_qval = score_plot_combined.qvals[i]
            elif mask > 0.0 and inblock:
                if max_logfc < mask:
                    max_logfc = mask
                    min_pval = score_plot_combined.pvals[i]
                    min_qval = score_plot_combined.qvals[i]
            elif mask <= 0.0 and inblock:
                inblock = False
                end = i
                blocks.append(Block(start + 1, end, end - start, max_logfc, 'increased_insertions', min_pval, min_qval))
                max_logfc = 0.0
                min_pval = 1.0
                min_qval = 1.0

        # Check for block at end
        if inblock:
            blocks.append(Block(start + 1, len(masking_plot_combined), len(masking_plot_combined) - start, max_logfc,
                                'increased_insertions', min_pval, min_qval))
        return blocks

    def decreased_insertions_blocks(self, masking_plot_combined, score_plot_combined):
        blocks = []
        inblock = False
        start = 0
        end = 0
        max_logfc = 0.0
        min_pval = 1.0
        min_qval = 1.0

        for i, mask in enumerate(masking_plot_combined):
            if mask < 0.0 and not inblock:
                inblock = True
                start = i
                max_logfc = mask
                min_pval = score_plot_combined.pvals[i]
                min_qval = score_plot_combined.qvals[i]
            elif mask < 0.0 and inblock:
                if max_logfc > mask:
                    max_logfc = mask
                    min_pval = score_plot_combined.pvals[i]
                    min_qval = score_plot_combined.qvals[i]
            elif mask >= 0.0 and inblock:
                inblock = False
                end = i
                blocks.append(Block(start + 1, end, end - start, max_logfc, 'decreased_insertions', min_pval, min_qval))
                max_logfc = 0.0
                min_pval = 1.0
                min_qval = 1.0

        # Check for block at end
        if inblock:
            blocks.append(Block(start + 1, len(masking_plot_combined), len(masking_plot_combined) - start, max_logfc,
                                'decreased_insertions', min_pval, min_qval))

        return blocks

    def peak_from_array(self, block_values):
        abs_max_value = max(numpy.absolute(block_values))
        if max(block_values) != abs_max_value:
            abs_max_value *= -1
        return abs_max_value

    # check if the reads are going in 1 particular direction, if so then use the max_logfc for that direction.
    def direction_for_block(self, block, forward_masking_plot, reverse_masking_plot):
        forward_max_logfc = self.peak_from_array(forward_masking_plot.combined[block.start - 1:block.end])
        reverse_max_logfc = self.peak_from_array(reverse_masking_plot.combined[block.start - 1:block.end])

        if numpy.absolute(forward_max_logfc) > numpy.absolute(reverse_max_logfc):
            if reverse_max_logfc == 0 or forward_max_logfc >= reverse_max_logfc + self.logfc_direction_change:
                block.max_logfc = forward_max_logfc

                return 'forward'
            else:
                return 'nodirection'
        elif numpy.absolute(forward_max_logfc) == numpy.absolute(reverse_max_logfc):
            return 'nodirection'
        else:
            if forward_max_logfc == 0 or reverse_max_logfc >= forward_max_logfc + self.logfc_direction_change:
                block.max_logfc = reverse_max_logfc
                return 'reverse'
            else:
                return 'nodirection'

    def merge_all_plots_choosing_peak_logfc(self, combined_plot, forward_plot, reverse_plot):
        genome_length = len(combined_plot.combined)
        masking_plot_combined = numpy.zeros(genome_length, dtype=int)
        for i in range(0, genome_length):
            peak_value_abs = max(
                numpy.absolute([combined_plot.combined[i], forward_plot.combined[i], reverse_plot.combined[i]]))
            if peak_value_abs == max([combined_plot.combined[i], forward_plot.combined[i], reverse_plot.combined[i]]):
                # its positive
                masking_plot_combined[i] = peak_value_abs
            else:
                # its negative
                masking_plot_combined[i] = peak_value_abs * -1

        return masking_plot_combined

    def block_generator(self):
        combined_plot = PlotParser(self.combined_mask_file)
        forward_masking_plot = PlotParser(self.forward_mask_file)
        reverse_masking_plot = PlotParser(self.reverse_mask_file)
        score_plot = ScoreParser(self.combined_score_file)

        masking_plot = self.merge_all_plots_choosing_peak_logfc(combined_plot, forward_masking_plot,
                                                                reverse_masking_plot)
        blocks = self.increased_insertions_blocks(masking_plot, score_plot) + self.decreased_insertions_blocks(
            masking_plot, score_plot)

        '''Filter out blocks which are less than the window size'''
        filtered_blocks = [b for b in blocks if b.block_length >= self.window_size]

        for b in filtered_blocks:
            b.direction = self.direction_for_block(b, forward_masking_plot, reverse_masking_plot)

        return filtered_blocks
