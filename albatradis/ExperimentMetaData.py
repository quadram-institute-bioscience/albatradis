
class ExperimentMetaData:
	def __init__(self, drug, target, detailed_target, impact, mic, induction, rep1_file, rep2_file):
		self.drug      = drug
		self.target    = target
		self.detailed_target = detailed_target
		self.impact    = impact
		self.mic       = mic
		self.induction = induction
		self.rep1_file = rep1_file
		self.rep2_file = rep2_file
