import subprocess
from tempfile import mkstemp
import os
import shutil
import csv

class TradisEssentiality:
	def __init__(self, tabfile, verbose, exec="tradis_essentiality.R", prefix="", plotnames="", analysis_type=""):
		self.tabfile = tabfile
		self.exec = exec
		self.verbose = verbose
		self.prefix = prefix
		self.plotnames = plotnames
		self.analysis_type = analysis_type

		fd, self.output_filename = mkstemp()
		fd, self.essential_filename = mkstemp()

	def construct_command(self):
		# print("Tabfile: "+ self.tabfile)
		return " ".join([self.exec,  self.tabfile])

	def run(self, plotname):
		cmd = self.construct_command()
		if self.verbose:
			print(cmd)
		subprocess.check_output(cmd, shell=True)
		
		self.replace_comma_tabs(self.tabfile +".all.csv", self.output_filename)
		shutil.copy(self.tabfile +".essen.csv", self.essential_filename)

		if self.verbose:
			print("all.csv\t" + self.output_filename)
			print("essen.csv\t" + self.essential_filename)

		if self.analysis_type == "original":
			condition_name = os.path.join(self.prefix, plotname + self.analysis_type + ".ess")
			shutil.copy(self.tabfile + ".essen.csv", condition_name)

		os.remove(self.tabfile +".all.csv")
		os.remove(self.tabfile +".essen.csv")
		os.remove(self.tabfile +".ambig.csv")
		os.remove(self.tabfile +".QC_and_changepoint_plots.pdf")
		os.remove(self.tabfile)
		
		return self
		
	def replace_comma_tabs(self, input_file, output_file):
		with open(input_file, newline='') as csvfile:
			comparison_reader = csv.reader(csvfile, delimiter=',')
			
			with open(output_file, 'w') as outputfh:
				for line in comparison_reader:
					outputfh.write("\t".join(line)+"\n")
