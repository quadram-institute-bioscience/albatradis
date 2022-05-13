'''Driver class'''
import logging
import os
import shutil
from tempfile import mkstemp

import quatradis.isp_analyse

from albatradis.BlockIdentifier import BlockIdentifier
from albatradis.Gene import Gene
from albatradis.GeneAnnotator import GeneAnnotator
from albatradis.PlotLog import PlotLog
from albatradis.PlotMasking import PlotMasking
from albatradis.PrepareEMBLFile import PrepareEMBLFile
from albatradis.PrepareInputFiles import PrepareInputFiles
from albatradis.TradisComparison import TradisComparison
from albatradis.TradisEssentiality import TradisEssentiality


class PlotEssentiality:
    def __init__(self, plotfile_obj, gene_insert_sites_filename, tradis_essentiality_filename, type,
                 only_essential_filename):
        self.plotfile_obj = plotfile_obj
        self.gene_insert_sites_filename = gene_insert_sites_filename
        self.tradis_essentiality_filename = tradis_essentiality_filename
        self.only_essential_filename = only_essential_filename
        self.type = type


class PlotAllEssentiality:
    def __init__(self, forward, reverse, combined, original, embl_filename):
        self.forward = forward
        self.reverse = reverse
        self.combined = combined
        self.original = original
        self.embl_filename = embl_filename


class BlockInsertions:
    def __init__(self, logger, plotfiles, plotnames, minimum_threshold, window_size, window_interval, verbose, minimum_logfc,
                 maximum_pvalue, maximum_qvalue, prefix, minimum_logcpm, minimum_block, span_gaps, emblfile,
                 report_decreased_insertions, strict_signal, use_annotation, prime_feature_size):
        self.logger = logger
        self.plotfiles = plotfiles
        self.plotnames = plotnames
        self.minimum_threshold = minimum_threshold
        self.window_size = window_size
        self.window_interval = window_interval
        self.verbose = verbose
        self.minimum_logfc = minimum_logfc
        self.maximum_pvalue = maximum_pvalue
        self.maximum_qvalue = maximum_qvalue
        self.prefix = prefix
        self.minimum_logcpm = minimum_logcpm
        self.minimum_block = minimum_block
        self.span_gaps = span_gaps
        self.emblfile = emblfile
        self.report_decreased_insertions = report_decreased_insertions
        self.strict_signal = strict_signal
        self.use_annotation = use_annotation
        self.prime_feature_size = prime_feature_size

        self.genome_length = 0
        self.annotation_file = ""
        self.output_plots = {}
        self.blocks = []
        self.genes = []

        if self.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.ERROR)

        if not os.path.exists(self.prefix):
            os.makedirs(self.prefix)

        if self.use_annotation:
            self.span_gaps = 0

    def run(self):
        plotfile_objects = self.prepare_input_files()
        essentiality_files = self.run_essentiality(plotfile_objects, self.plotnames)

        forward_plotfile, reverse_plotfile, combined_plotfile, original_plotfile,\
            forward_scoresfile, reverse_scoresfile, combined_scoresfile, original_scoresfile = self.run_comparisons(essentiality_files)
        self.output_plots = self.mask_plots(combined_plotfile)
        self.genes = self.gene_statistics(forward_plotfile, reverse_plotfile, combined_plotfile,
                                          combined_scoresfile, self.window_size)
        self.cleanup(plotfile_objects, essentiality_files)

        return self

    def prepare_input_files(self):
        plotfile_objects = {}

        print("1. Preparing input files.\n")
        if self.verbose:
            print("1.1. Preparing EMBL file.\n")

        self.annotation_file = PrepareEMBLFile(self.plotfiles[0], self.minimum_threshold, self.window_size,
                                               self.window_interval, self.use_annotation, self.prime_feature_size,
                                               self.emblfile).create_file()

        if self.verbose:
            print(" Embl:\t" + self.annotation_file + "\n")
            print("1.2. Split input plot files into forward, reverse, combined and original.\n")

        i = 1
        n = len(self.plotfiles)

        for plotfile in self.plotfiles:
            p = PrepareInputFiles(plotfile, self.minimum_threshold)
            p.create_all_files()
            p.embl_filename = self.annotation_file
            plotfile_objects[plotfile] = p

            if self.verbose:
                print("Split input plot " + plotfile + ".\n")
                print("Forward plot:\t" + p.forward_plot_filename)
                print("reverse plot:\t" + p.reverse_plot_filename)
                print("combined plot:\t" + p.combined_plot_filename)
                print("original plot:\t" + p.original_plot_filename)

            print("\t" + str(i) + " of " + str(n) + " files processed.\n")
            i = i + 1
            self.genome_length = p.genome_length()
        return plotfile_objects

    def analyse_insert_sites(self, embl_filename, plotfile):
        fd, output_filename = mkstemp()
        quatradis.isp_analyse.analyse_insert_sites(embl_filename, [plotfile])
        plotfile_prefix = os.path.basename(plotfile)
        plotfile_prefix = plotfile_prefix.split(sep=".")[0]
        shutil.copy(plotfile_prefix + ".tradis_gene_insert_sites.csv", output_filename)
        os.remove(plotfile_prefix + ".tradis_gene_insert_sites.csv")
        return output_filename

    def essentiality(self, embl_filename, plotfile_objects, plotfile, plotname, filetype):

        analysis_file = self.analyse_insert_sites(embl_filename, getattr(plotfile_objects[plotfile], filetype + "_plot_filename"))

        e = TradisEssentiality(analysis_file, self.verbose, prefix=self.prefix, plotnames=self.plotnames, analysis_type=filetype)
        e.run(plotname)
        pe = PlotEssentiality(plotfile, analysis_file, e.output_filename, filetype, e.essential_filename)

        if self.verbose:
            print("Essentiality:\t" + filetype + "\t" + e.output_filename)

        return pe

    def run_essentiality(self, plotfile_objects, plotnames):
        essentiality_files = {}

        n = 4 * len(plotfile_objects)
        i = 4
        j = 0

        print("2. Find essential genes for each split file.\n")

        for plotfile in plotfile_objects:
            f = self.essentiality(plotfile_objects[plotfile].embl_filename, plotfile_objects, plotfile, plotnames[j], 'forward')
            r = self.essentiality(plotfile_objects[plotfile].embl_filename, plotfile_objects, plotfile, plotnames[j], 'reverse')
            c = self.essentiality(plotfile_objects[plotfile].embl_filename, plotfile_objects, plotfile, plotnames[j], 'combined')
            o = self.essentiality(self.emblfile, plotfile_objects, plotfile, plotnames[j], 'original')
            k = plotfile_objects[plotfile].embl_filename
            essentiality_files[plotfile] = PlotAllEssentiality(f, r, c, o, k)
            print("\t" + str(i) + " of " + str(n) + " files processed.\n")
            i = i + 4
            j = j + 1

        return essentiality_files

    def run_comparisons(self, essentiality_files):

        print("3. Run comparison between condition and control and write plot files and gene report.\n")

        forward_plotfile, forward_scoresfile = self.generate_logfc_plot('forward', essentiality_files)
        print("Forward insertions have been compared.\n")
        reverse_plotfile, reverse_scoresfile = self.generate_logfc_plot('reverse', essentiality_files)
        print("Reverse insertions have been compared.\n")
        combined_plotfile, combined_scoresfile = self.generate_logfc_plot('combined', essentiality_files)
        print("Combined insertions have been compared.\n")
        original_plotfile, original_scoresfile = self.generate_logfc_plot('original', essentiality_files)
        print("Original insertions have been compared.\n")

        return forward_plotfile, reverse_plotfile, combined_plotfile, original_plotfile, \
               forward_scoresfile, reverse_scoresfile, combined_scoresfile, original_scoresfile

    def generate_logfc_plot(self, analysis_type, essentiality_files):
        files = [getattr(essentiality_files[plotfile], analysis_type).tradis_essentiality_filename for plotfile in
                 self.plotfiles]

        only_ess_files = [getattr(essentiality_files[plotfile], analysis_type).only_essential_filename for plotfile in
                          self.plotfiles]

        annotation_files = [essentiality_files[plotfile].embl_filename for plotfile in self.plotfiles]
        mid = int(len(files) / 2)

        t = TradisComparison(files[:mid], files[mid:], self.verbose, self.minimum_block, only_ess_files[:mid],
                             only_ess_files[mid:], analysis_type, self.prefix)
        t.run()
        p = PlotLog(t.output_filename, self.genome_length, self.minimum_logfc, self.maximum_pvalue, self.maximum_qvalue,
                    self.minimum_logcpm, self.window_size, self.span_gaps, self.report_decreased_insertions,
                    annotation_files[0])
        p.construct_plot_file()

        # Move temp files to output dir
        renamed_csv_file = os.path.join(self.prefix, analysis_type + ".csv")
        renamed_plot_file = os.path.join(self.prefix, analysis_type + ".logfc.plot")
        renamed_scores_file = os.path.join(self.prefix, analysis_type + ".pqvals.plot")
        shutil.move(t.output_filename, renamed_csv_file)
        shutil.move(p.plot_filename, renamed_plot_file)
        shutil.move(p.scores_filename, renamed_scores_file)

        if self.verbose:
            print("Comparison:\t" + renamed_csv_file)
            print("Plot log:\t" + renamed_plot_file)
        return renamed_plot_file, renamed_scores_file

    def merge_windows(self, windows):

        start_window = windows[0]
        i = 1
        merged_windows = []
        while i < len(windows):
            next_window = windows[i]

            if (next_window.categories[0] == "increased_mutants_at_end_of_gene" or \
                    next_window.categories[0] == "increased_mutants_at_start_of_gene" or \
                    next_window.categories[0] == "decreased_mutants_at_end_of_gene" or \
                    next_window.categories[0] == "decreased_mutants_at_start_of_gene"):
                next_window.categories[0] = "knockout"
            if (start_window.categories[0] == "increased_mutants_at_end_of_gene" \
                    or start_window.categories[0] == "increased_mutants_at_start_of_gene" or \
                    start_window.categories[0] == "decreased_mutants_at_end_of_gene" or \
                    start_window.categories[0] == "decreased_mutants_at_start_of_gene"):
                start_window.categories[0] = "knockout"
            category_equal = next_window.categories[0] == start_window.categories[0]
            expression_equal = next_window.expression_from_blocks() == start_window.expression_from_blocks()
            direction_equal = next_window.direction_from_blocks() == start_window.direction_from_blocks()
            if next_window.start <= start_window.end and category_equal and expression_equal and direction_equal:
                start_window.end = next_window.end
                start_window.max_logfc = max(start_window.max_logfc, next_window.max_logfc)
                start_window.min_pvalue = min(start_window.min_pvalue, next_window.min_pvalue)
                start_window.min_qvalue = min(start_window.min_qvalue, next_window.min_qvalue)
            else:
                copy = Gene(start_window.feature, start_window.blocks)
                copy.start = start_window.start
                copy.end = start_window.end
                copy.max_logfc = start_window.max_logfc
                copy.min_pvalue = start_window.min_pvalue
                copy.min_qvalue = start_window.min_qvalue
                copy.gene_name = str(copy.start) + "_" + str(copy.end)
                copy.categories.append(start_window.categories[0])
                merged_windows.append(copy)
                start_window = next_window

            i += 1

        merged_windows.append(start_window)
        return merged_windows

    def gene_statistics(self, forward_plotfile, reverse_plotfile, combined_plotfile, combined_scorefile, window_size):
        b = BlockIdentifier(combined_plotfile, forward_plotfile, reverse_plotfile, combined_scorefile, window_size)
        blocks = b.block_generator()
        annotationfile = self.emblfile
        if self.use_annotation:
            annotationfile = self.annotation_file

        genes = GeneAnnotator(self.annotation_file, blocks).annotate_genes()
        intergenic_blocks = [block for block in blocks if block.intergenic]

        if not self.use_annotation:
            all_genes = self.merge_windows(genes)
        else:
            all_genes = []
            for g in genes:
                all_genes.append(g)

        if len(all_genes) == 0 and len(intergenic_blocks) == 0:
            print("No significant genes found for chosen parameters.\n")
            return []

        self.write_gene_report(all_genes, intergenic_blocks)
        self.write_regulated_gene_report(all_genes, intergenic_blocks)

        # if self.verbose:
        # self.print_genes_intergenic(genes,intergenic_blocks)

        return genes

    def write_gene_report(self, genes, intergenic_blocks):
        block_filename = os.path.join(self.prefix, "gene_report.csv")

        with open(block_filename, 'w') as bf:
            bf.write(str(genes[0].header()) + "\n")
            if not self.use_annotation:
                for i in genes:
                    bf.write(i.window_string() + "\n")
            else:
                for i in genes:
                    bf.write(str(i) + "\n")
            for b in intergenic_blocks:
                bf.write(str(b) + "\n")

    def write_regulated_gene_report(self, genes, intergenic_blocks):
        regulated_genes = [g for g in genes if g.category() == 'upregulated' or g.category() == 'downregulated']
        if len(regulated_genes) > 0:
            block_filename = os.path.join(self.prefix, "regulated_gene_report.csv")
            with open(block_filename, 'w') as bf:
                bf.write(str(regulated_genes[0].header()) + "\n")
                for i in regulated_genes:
                    bf.write(str(i) + "\n")

    def print_genes_intergenic_blocks(self, genes, intergenic_blocks):
        print(genes[0].header())
        for i in genes:
            print(i)
        for b in intergenic_blocks:
            print(b)

    def mask_plots(self, combined_plotfile):
        pm = PlotMasking(self.plotfiles, combined_plotfile, self.strict_signal)
        renamed_plot_files = {}

        for pfile in pm.output_plot_files:
            original_basefile = os.path.join(self.prefix, os.path.basename(pfile))
            renamed_file = original_basefile.replace('.gz', '')
            shutil.copy(pm.output_plot_files[pfile], renamed_file)
            os.remove(pm.output_plot_files[pfile])
            renamed_plot_files[pfile] = renamed_file

            if self.verbose:
                print("Masked: " + renamed_file)
        return renamed_plot_files

    def cleanup(self, plotfile_objects, essentiality_files):

        # initial plot files
        for p in plotfile_objects.values():
            os.remove(p.forward_plot_filename)
            os.remove(p.reverse_plot_filename)
            os.remove(p.combined_plot_filename)
            # os.remove(p.original_plot_filename)
            if os.path.exists(p.embl_filename):
                shutil.move(p.embl_filename, os.path.join(self.prefix, "annotation.embl"))
