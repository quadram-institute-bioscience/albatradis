from Bio import SeqIO

class EMBLReader:
	def __init__(self, filename):
		self.filename = filename
		self.features_to_ignore = ['source','gene']

	def read_annotation_features(self):
		record =  SeqIO.read(self.filename, "embl")
		return [f for f in record.features if f.type not in self.features_to_ignore]
	
