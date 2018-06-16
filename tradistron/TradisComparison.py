import subprocess
from tempfile import mkstemp
import os

class TradisComparison:
	def __init__(self, condition_files, control_files, verbose, minimum_block, exec="tradis_comparison.R"):
		self.condition_files = condition_files
		self.control_files = control_files
		self.exec     = exec
		self.verbose  = verbose
		self.minimum_block = minimum_block
		
		fd, self.output_filename = mkstemp()
		fd, self.conditions_fofn = mkstemp()
		fd, self.controls_fofn = mkstemp()

	def create_fofn(self):
		with open(self.conditions_fofn, 'w') as fileh:
			if len(self.condition_files) == 1:
				fileh.write(self.condition_files[0]+"\n")
				fileh.write(self.condition_files[0]+"\n")
			else:
				for i in self.condition_files:
					fileh.write(i + "\n")
					
		with open(self.controls_fofn, 'w') as fileh:
			if len(self.control_files) == 1:
				fileh.write(self.control_files[0]+"\n")
				fileh.write(self.control_files[0]+"\n")
			else:
				for i in self.control_files:
					fileh.write(i + "\n")
		return self
		
	def construct_command(self):
		return " ".join([self.exec, '-f', '-t', str(self.minimum_block), '--controls', self.controls_fofn, '--conditions', self.conditions_fofn])
		
	def run(self):
		self.create_fofn()
		
		if self.verbose:
			print(self.construct_command())
		subprocess.check_output(self.construct_command(), shell=True)
		os.rename(".output.csv", self.output_filename)

		return self
	
