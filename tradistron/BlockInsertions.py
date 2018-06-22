'''Driver class'''
from tradistron.PrepareInputFiles import PrepareInputFiles
import logging
import os
import sys
import time
from tradistron.TradisGeneInsertSites import TradisGeneInsertSites
from tradistron.PrepareInputFiles     import PrepareInputFiles
from tradistron.TradisEssentiality    import TradisEssentiality
from tradistron.TradisComparison      import TradisComparison
from tradistron.PlotLog               import PlotLog
from tradistron.PlotMasking           import PlotMasking
from tradistron.BlockIdentifier       import BlockIdentifier
from tradistron.GeneAnnotator        import GeneAnnotator

class PlotEssentiality:
	def __init__(self, plotfile_obj,gene_insert_sites_filename, tradis_essentiality_filename, type, only_essential_filename):
		self.plotfile_obj = plotfile_obj
		self.gene_insert_sites_filename = gene_insert_sites_filename
		self.tradis_essentiality_filename = tradis_essentiality_filename
		self.only_essential_filename = only_essential_filename
		self.type = type
		
class PlotAllEssentiality:
	def __init__(self, forward, reverse, combined):
		self.forward = forward
		self.reverse = reverse
		self.combined = combined

class BlockInsertions:
	def __init__(self, logger,plotfiles, minimum_threshold, window_size, window_interval, verbose, minimum_logfc, pvalue, prefix, minimum_logcpm, minimum_block,span_gaps, emblfile):
		self.logger            = logger
		self.plotfiles         = plotfiles
		self.minimum_threshold = minimum_threshold
		self.window_size       = window_size
		self.window_interval   = window_interval
		self.verbose           = verbose
		self.minimum_logfc     = minimum_logfc
		self.pvalue            = pvalue
		self.prefix            = prefix
		self.minimum_logcpm    = minimum_logcpm
		self.minimum_block     = minimum_block
		self.span_gaps         = span_gaps
		self.emblfile          = emblfile  
		
		self.genome_length = 0
		self.forward_plotfile = ""
		self.reverse_plotfile = ""
		self.combined_plotfile = ""
		self.output_plots = {}
		self.blocks = []
		
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
			
		if not os.path.exists(self.prefix ):
			os.makedirs(self.prefix )
		
	def run(self):
		plotfile_objects = self.prepare_input_files()
		essentiality_files = self.run_essentiality(plotfile_objects)
		
		self.run_comparisons(essentiality_files)
		self.output_plots = self.mask_plots()
		self.genes = self.gene_statistics(self.forward_plotfile, self.reverse_plotfile, self.combined_plotfile, self.window_size)
		
		return self
		
	def prepare_input_files(self):
		plotfile_objects = {}
		for plotfile in self.plotfiles:
			p = PrepareInputFiles(plotfile, self.minimum_threshold, self.window_size, self.window_interval )
			p.create_all_files()
			plotfile_objects[plotfile] = p
			
			if self.verbose:
				print("Forward plot:\t" + p.forward_plot_filename)
				print("reverse plot:\t" + p.reverse_plot_filename)
				print("combined plot:\t" + p.combined_plot_filename)
				print("Embl:\t" + p.embl_filename)
			
			self.genome_length = p.genome_length()
		return plotfile_objects
	
	def essentiality(self, plotfile_objects, plotfile, filetype):
		g = TradisGeneInsertSites(plotfile_objects[plotfile].embl_filename, getattr(plotfile_objects[plotfile], filetype + "_plot_filename"), self.verbose)
		g.run()
		e = TradisEssentiality(g.output_filename, self.verbose)
		e.run()
		pe = PlotEssentiality(plotfile, g.output_filename, e.output_filename, filetype, e.essential_filename)
		
		if self.verbose:
			print("Essentiality:\t" + filetype + "\t" + e.output_filename)
		return pe
		
	def run_essentiality(self, plotfile_objects):
		essentiality_files = {}
		for plotfile in plotfile_objects:
			f = self.essentiality(plotfile_objects, plotfile, 'forward')
			r = self.essentiality(plotfile_objects, plotfile, 'reverse')
			c = self.essentiality(plotfile_objects, plotfile, 'combined')
			essentiality_files[plotfile] = PlotAllEssentiality(f,r,c)

		return essentiality_files
		
	def run_comparisons(self, essentiality_files):
		self.forward_plotfile = self.generate_logfc_plot('forward',essentiality_files)
		self.reverse_plotfile = self.generate_logfc_plot('reverse',essentiality_files)
		self.combined_plotfile = self.generate_logfc_plot('combined',essentiality_files)
			
	def generate_logfc_plot(self, analysis_type, essentiality_files):
		files = [getattr(essentiality_files[plotfile], analysis_type).tradis_essentiality_filename for plotfile in self.plotfiles]
		
		only_ess_files = [getattr(essentiality_files[plotfile], analysis_type).only_essential_filename for plotfile in self.plotfiles]
		
		mid = int(len(files)  / 2)
		
		t = TradisComparison(files[:mid],files[mid:], self.verbose, self.minimum_block, only_ess_files[:mid], only_ess_files[mid:])
		t.run()
		p = PlotLog(t.output_filename, self.genome_length, self.minimum_logfc, self.pvalue, self.minimum_logcpm, self.window_size, self.span_gaps)
		p.construct_plot_file()
		renamed_csv_file  = os.path.join(self.prefix, analysis_type + ".csv")
		renamed_plot_file = os.path.join(self.prefix, analysis_type + ".plot")
		
		os.rename(t.output_filename, renamed_csv_file)
		os.rename(p.output_filename, renamed_plot_file)
		if self.verbose:
			print("Comprison:\t"+ renamed_csv_file)
			print("Plot log:\t"+ renamed_plot_file)
		return renamed_plot_file
		
		
	def gene_statistics(self,forward_plotfile, reverse_plotfile, combined_plotfile, window_size):
		b = BlockIdentifier(combined_plotfile, forward_plotfile, reverse_plotfile, window_size)
		blocks = b.block_generator()
		genes = GeneAnnotator(self.emblfile, blocks).annotate_genes()
		intergenic_blocks = [block for block in blocks if block.intergenic]
		
		if len(genes) == 0:
			return []
		
		block_filename = os.path.join(self.prefix, "gene_report.csv")
		with open(block_filename, 'w') as bf:
			bf.write(str(genes[0].header())+"\n")
			for i in genes:
				bf.write(str(i)+"\n")
				
			for b in intergenic_blocks:
				bf.write(str(b)+"\n")
				
		if self.verbose:
			print(genes[0].header())		
			for i in genes:
				print(i)
				
			for b in intergenic_blocks:
				print(b)
		
		return genes
		
	def mask_plots(self):
		pm = PlotMasking(self.plotfiles, self.combined_plotfile )
		renamed_plot_files = {}
		
		for pfile in pm.output_plot_files:
			original_basefile  = os.path.join(self.prefix, os.path.basename(pfile) )
			renamed_file = original_basefile.replace('.gz','')
			os.rename(pm.output_plot_files[pfile], renamed_file)
			renamed_plot_files[pfile] = renamed_file
			
			if self.verbose:
				print("Masked: " + renamed_file )
		return renamed_plot_files
		
		
		
		
	
		