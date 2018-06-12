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

class PlotEssentiality:
	def __init__(self, plotfile_obj,gene_insert_sites_filename, tradis_essentiality_filename, type):
		self.plotfile_obj = plotfile_obj
		self.gene_insert_sites_filename = gene_insert_sites_filename
		self.tradis_essentiality_filename = tradis_essentiality_filename
		self.type = type
		
class PlotAllEssentiality:
	def __init__(self, forward, reverse, combined):
		self.forward = forward
		self.reverse = reverse
		self.combined = combined

class BlockInsertions:
	def __init__(self, logger,plotfiles, minimum_threshold, window_size, window_interval, verbose, minimum_logfc, pvalue, prefix, minimum_logcpm):
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
		
		self.genome_length = 0
		self.combined_plotfile = ""
		self.output_plots = {}
		
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
		pe = PlotEssentiality(plotfile, g.output_filename, e.output_filename, filetype)
		
		if self.verbose:
			print("Essentiality:\t" + filetype + "\t" + e.output_filename)
		return pe
		
	def run_essentiality(self,plotfile_objects):
		essentiality_files = {}
		for plotfile in plotfile_objects:
			f = self.essentiality(plotfile_objects, plotfile, 'forward')
			r = self.essentiality(plotfile_objects, plotfile, 'reverse')
			c = self.essentiality(plotfile_objects, plotfile, 'combined')
			essentiality_files[plotfile] = PlotAllEssentiality(f,r,c)

		return essentiality_files
		
	def run_comparisons(self, essentiality_files):
		files = [essentiality_files[plotfile].forward.tradis_essentiality_filename for plotfile in essentiality_files]
		t = TradisComparison([files[0]],[files[1]], self.verbose)
		t.run()
		p = PlotLog(t.output_filename, self.genome_length, self.minimum_logfc, self.pvalue, self.minimum_logcpm)
		p.construct_plot_file()
		os.rename(t.output_filename, os.path.join(self.prefix,"forward.csv"))
		os.rename(p.output_filename, os.path.join(self.prefix,"forward.plot"))
		if self.verbose:
			print("Comprison\t"+ os.path.join(self.prefix,"forward.csv"))
			print("Plot log:\t"+ os.path.join(self.prefix,"forward.plot"))
		
		files = [essentiality_files[plotfile].reverse.tradis_essentiality_filename for plotfile in essentiality_files]
		t = TradisComparison([files[0]],[files[1]], self.verbose)
		t.run()
		p = PlotLog(t.output_filename, self.genome_length, self.minimum_logfc, self.pvalue, self.minimum_logcpm)
		p.construct_plot_file()
		os.rename(t.output_filename, os.path.join(self.prefix,"reverse.csv"))
		os.rename(p.output_filename, os.path.join(self.prefix,"reverse.plot"))
		if self.verbose:
			print("Comprison\t"+ os.path.join(self.prefix,"reverse.csv"))
			print("Plot log:\t"+ os.path.join(self.prefix,"reverse.plot"))
		
		files = [essentiality_files[plotfile].combined.tradis_essentiality_filename for plotfile in essentiality_files]
		t = TradisComparison([files[0]],[files[1]], self.verbose)
		t.run()
		p = PlotLog(t.output_filename, self.genome_length, self.minimum_logfc, self.pvalue, self.minimum_logcpm)
		p.construct_plot_file()
		os.rename(t.output_filename, os.path.join(self.prefix,"combined.csv"))
		os.rename(p.output_filename, os.path.join(self.prefix,"combined.plot"))
		self.combined_plotfile = os.path.join(self.prefix,"combined.plot")
		if self.verbose:
			print("Comprison\t"+ os.path.join(self.prefix,"combined.csv"))
			print("Plot log:\t"+ os.path.join(self.prefix,"combined.plot"))
		
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
		
		
		
		
	
		