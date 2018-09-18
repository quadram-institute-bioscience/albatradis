import re
import os
import numpy

class ExperimentCollection:
	def __init__(self, experiments, controls, reference):
		self.experiments = experiments
		self.controls  = controls
		self.reference = reference
		# 'target', 'detailed_target', 'impact',
		self.property_names = ['drug',  'mic', 'induction']
		self.in_common = self.properties_in_common()
		self.project_name = self.generate_project_name(self.in_common)
		
	def properties_in_common(self):
		properties = {}
		for p in self.property_names:
			property_to_values = numpy.unique([getattr(e, p) for e in self.experiments ])
			if len(property_to_values) > 1:
				continue
			else:
				properties[p] = property_to_values[0]
		return properties
		
	def generate_project_name(self, properties):
		name_parts = []
		for p in self.property_names:
			if p in properties:
				name_parts.append(str(properties[p]))
		
		merged_name = "_".join(name_parts)
		regex = re.compile('[\W]')
		cleaned_name = regex.sub('_', merged_name)
		return cleaned_name
		
	def properties_not_in_common(self):
		return [ x for x in self.property_names if x not in self.in_common ]

	def sorted_experiments(self):
		variable_attributes  = self.properties_not_in_common()
		
		# TODO: make this work more generally
		if 'mic' in variable_attributes and 'induction' in variable_attributes:
			return sorted(self.experiments, key = lambda x: (x.mic, x.induction))
		elif 'mic' in variable_attributes:
			return sorted(self.experiments, key = lambda x: (x.mic))
		elif 'induction' in variable_attributes:
			return sorted(self.experiments, key = lambda x: (x.induction))
		else:
			return sorted(self.experiments, key = lambda x:(x.rep1_file))

	def __str__(self):
		output  = "project." + self.project_name + ".sequence=" + self.reference + "\n"
		output += "project." + self.project_name + ".userplot=" + " ".join(self.controls)
		
		for e in self.sorted_experiments():
			if e.rep1_file == 'missing' or e.rep2_file == 'missing':
				continue
			if e.rep1_file != None and os.path.exists(e.rep1_file):
				output += " " + e.rep1_file
			if e.rep2_file != None and os.path.exists(e.rep2_file):
		 		output += " " + e.rep2_file
		
		output += "\n"
		return output
		