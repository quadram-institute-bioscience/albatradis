import subprocess
from tempfile import mkstemp
import os

class TradisEssentiality:
	def __init__(self, tabfile, verbose, exec="tradis_essentiality.R"):
		self.tabfile = tabfile
		self.exec     = exec
		self.verbose = verbose
		
		fd, self.output_filename = mkstemp()
		fd, self.essential_filename = mkstemp()
		fd, self.ambig_filename = mkstemp()

	def construct_command(self):
		return " ".join([self.exec,  self.tabfile])
		
	def run(self):
		if self.verbose:
			print(self.construct_command())
		subprocess.check_output(self.construct_command(), shell=True)
		os.rename(self.tabfile +".all.csv", self.output_filename)
		os.rename(self.tabfile +".essen.csv", self.essential_filename)
		os.rename(self.tabfile +".ambig.csv", self.ambig_filename)
		
		return self
