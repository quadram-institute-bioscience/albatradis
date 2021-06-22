import subprocess
from tempfile import mkstemp
import os
import shutil

class TradisGeneInsertSites:
	def __init__(self, emblfile, plotfile, verbose, exec="tradis_gene_insert_sites"):
		self.emblfile = emblfile
		self.plotfile = plotfile
		self.exec     = exec
		self.verbose  = verbose
		
		fd, self.output_filename = mkstemp()

	def construct_command(self):
		return " ".join([self.exec,  self.emblfile, self.plotfile ])
		
	def run(self):
		cmd = self.construct_command()
		if self.verbose:
			print(cmd)
		subprocess.check_output(cmd, shell=True)
		plotfile_prefix = os.path.basename(self.plotfile)
		plotfile_prefix = plotfile_prefix.split(sep=".")[0]
		shutil.copy(plotfile_prefix +".tradis_gene_insert_sites.csv", self.output_filename)
		os.remove(plotfile_prefix +".tradis_gene_insert_sites.csv")
		return self
		
