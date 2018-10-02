from Bio import SeqIO

class EMBLReader:
	def __init__(self, filename):
		self.filename = filename
		self.features_to_ignore = ['source','gene']
		self.genome_length = 0

	def read_annotation_features(self):
		record =  SeqIO.read(self.filename, "embl")
		self.genome_length = len(record.seq)
		
		return [f for f in record.features if f.type not in self.features_to_ignore]
	
