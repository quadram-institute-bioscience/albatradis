from Bio import SeqIO
from shutil import copyfile

class EMBLReader:
	def __init__(self, filename):
		self.filename = filename
		self.features_to_ignore = ['source', 'gene']
		self.genome_length = 0
		self.features = self.read_annotation_features()
		self.genes_to_features = self.gene_names_to_features()

	def read_annotation_features(self):
		self.record =  SeqIO.read(self.filename, "embl")
		self.genome_length = len(self.record.seq)

		return [f for f in self.record.features if f.type not in self.features_to_ignore]
	
	def gene_names_to_features(self):
		genes_to_features = {}
		for f in self.features:
			gene_name = self.feature_to_gene_name(f)
			genes_to_features[gene_name] = f
		return genes_to_features
			
	def feature_to_gene_name(self, feature):
		gene_name_val = str(feature.location.start) + "_" + str(feature.location.end)
		if "gene" in feature.qualifiers:
			gene_name_val = feature.qualifiers["gene"][0]
		return gene_name_val

