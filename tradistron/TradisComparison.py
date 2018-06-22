import subprocess
from tempfile import mkstemp
import os

class GeneEssentiality:
	def __init__(self):
		self.condition = 0
		self.control = 0
		
	def status(self):
		if self.condition > self.control:
			return 'changed_to_essential'
		elif self.condition < self.control:
			return 'changed_to_nonessential'
		elif self.condition == self.control and self.condition > 0:
			return 'always_essential'
		elif self.condition == 0 and self.control == 0:
			return 'always_nonessential'
		else:
			return 'unknown'

class TradisComparison:
	def __init__(self, condition_files, control_files, verbose, minimum_block, only_ess_files_condition, only_ess_files_control, exec="tradis_comparison.R"):
		self.condition_files = condition_files
		self.control_files = control_files
		self.exec     = exec
		self.verbose  = verbose
		self.minimum_block = minimum_block
		self.only_ess_files_condition = only_ess_files_condition
		self.only_ess_files_control = only_ess_files_control
		
		fd, self.output_filename = mkstemp()
		fd, self.conditions_fofn = mkstemp()
		fd, self.controls_fofn = mkstemp()

	
	# add an extra column to the end to indicate if a region has changed from
	# essential to non-essential, or vica versa.
	# if there is discordance between the cases or controls, mark as ambiguous

	def gene_names_from_essentiality_file(self,filename):
		gene_names = []
		with open(filename, 'r') as fileh:
			rows = [line.split(',') for line in fileh]
			
			gene_names = [ r[1] for r in rows if r[1] != 'gene_name']
		return gene_names
		
	def all_gene_essentiality(self):
		all_gene_names = self.gene_names_from_essentiality_file(self.output_filename)
		
		genes_ess = {g: GeneEssentiality() for g in all_gene_names}

		for f in self.only_ess_files_condition:
			ess_gene_names = self.gene_names_from_essentiality_file(f)
			for e in ess_gene_names:
				genes_ess[e].condition += 1
				
		for f in self.only_ess_files_control:
			ess_gene_names = self.gene_names_from_essentiality_file(f)
			for e in ess_gene_names:
				genes_ess[e].control += 1
		return genes_ess
	
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
		
	def add_gene_essentiality_to_file(self, input_filename,output_filename, genes_ess):
		with open(input_filename, 'r') as inputfh:
			output_content = []
			input_content = inputfh.read().splitlines()
			for i,line in enumerate(input_content):
				cell = line.split("\t")
				if i == 0:
					cell.append("Essentiality")
				elif cell[0] in genes_ess:
					cell.append(genes_ess[cell[0]].status)
				else:
					cell.append('unknown')
				output_content.append(cell)
					
			with open(output_filename, 'w') as outputfh:
				for line in output_content:
					outputfh.write("\t".join(line)+"\n")
		return self

		
	def construct_command(self):
		return " ".join([self.exec, '-f', '-t', str(self.minimum_block), '--controls', self.controls_fofn, '--conditions', self.conditions_fofn])
		
	def run(self):
		self.create_fofn()
		
		if self.verbose:
			print(self.construct_command())
		subprocess.check_output(self.construct_command(), shell=True)
		genes_ess = self.all_gene_essentiality()
		self.add_gene_essentiality_to_file(".output.csv", self.output_filename, genes_ess)
		

		return self
	
