import os
import sys
import re
import csv
import numpy

from tradistron.PresenceAbsence import PresenceAbsence
from tradistron.ExperimentMetaData import ExperimentMetaData
from tradistron.ExperimentCollection import ExperimentCollection

class ArtemisProject:
	def __init__(self, outputfile, verbose, experiments_metadata, reference, controls):
		self.outputfile           = outputfile
		self.verbose              = verbose
		self.experiments_metadata = experiments_metadata
		self.reference            = reference
		self.controls             = controls
		
		self.experiments = self.populate_experiment_meta_data()

	def populate_experiment_meta_data(self):
		experiments = []
		with open(self.experiments_metadata, 'r') as fileh:
			reader = csv.reader(fileh, delimiter=",")
			for r in reader:
				if len(r) != 8:
					continue
				if r[0] == 'Drug':
					continue
				experiments.append(ExperimentMetaData(r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7]))
		return experiments
		
	# takes in an array of Experiments and an attribute to split them on
	def split_experiments_by_attr(self, experiments, split_attr):
		experiments_by_attr = {}
		split_vals = numpy.unique([getattr(e, split_attr) for e in experiments])
		for d in split_vals:
			experiments_by_attr[d] = [e for e in experiments if getattr(e, split_attr) == d]
		return experiments_by_attr
	
	def create_project_file(self):
		# ['drug', 'target', 'detailed_target', 'impact', 'mic', 'induction']
		with open(self.outputfile, 'w') as fh:
			
			# 1 Drug, all MIC and induction
			drug_exp = self.split_experiments_by_attr(self.experiments, 'drug')
			for f in drug_exp:
				ec_drug = ExperimentCollection(drug_exp[f], self.controls, self.reference)
				fh.write(str(ec_drug))
				
				# 1 Drug 1 MIC
				mic_drug_exp = self.split_experiments_by_attr(drug_exp[f], 'mic')
				for m in mic_drug_exp:
					ec = ExperimentCollection(mic_drug_exp[m], self.controls, self.reference)
					fh.write(str(ec))

				# 1 Drug 1 induction
				ind_drug_exp = self.split_experiments_by_attr(drug_exp[f], 'induction')
				for m in ind_drug_exp:
					ec = ExperimentCollection(ind_drug_exp[m], self.controls, self.reference)
					fh.write(str(ec))
				
			# all Drug, 1 MIC, all induction
			mic_exp = self.split_experiments_by_attr(self.experiments, 'mic')
			for f in mic_exp:
				ec = ExperimentCollection(mic_exp[f], self.controls, self.reference)
				fh.write(str(ec))
				
				# 1 MIC 1 induction
				ind_drug_exp = self.split_experiments_by_attr(mic_exp[f], 'induction')
				for m in ind_drug_exp:
					ec = ExperimentCollection(ind_drug_exp[m], self.controls, self.reference)
					fh.write(str(ec))
			
			# all Drug, all MIC, 1 induction
			induction_exp = self.split_experiments_by_attr(self.experiments, 'induction')
			for f in induction_exp:
				ec = ExperimentCollection(induction_exp[f], self.controls, self.reference)
				fh.write(str(ec))

		return self
