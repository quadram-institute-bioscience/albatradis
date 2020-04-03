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
		if self.verbose:
			print(self.construct_command())
		subprocess.check_output(self.construct_command(), shell=True)
		file = os.path.basename(self.plotfile)
		file1 = file.split(sep =  ".")[0]
		print(file1 +".tradis_gene_insert_sites.csv")
		print(os.getcwd())
		shutil.copy(file1 +".tradis_gene_insert_sites.csv", self.output_filename)
		os.remove(file1 +".tradis_gene_insert_sites.csv")
		return self
		
