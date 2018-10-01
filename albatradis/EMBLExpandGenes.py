''' Given an annotation file, take each gene, and create a new feature at the start and end to capture promotors'''
from albatradis.EMBLReader import EMBLReader

class FeatureProperties:
	def __init__(self, 	start, end, direction, gene_name):
		self.start = start
		self.end = end
		self.direction = direction
		self.gene_name = gene_name

class EMBLExpandGenes:
	def __init__(self, embl_file, feature_size):
		self.embl_file = embl_file
		self.feature_size = feature_size
		self.features = EMBLReader(self.embl_file).read_annotation_features()
	
	def create_3_5_prime_features(self):
		new_features = []
		for feature in self.features:
			gene_name = self.feature_to_gene_name(feature)
			
			# The gene itself
			new_features.append(FeatureProperties(f.location.start, f.location.end, feature.strand, gene_name))
			
			# forward direction
			if feature.strand == 1:
				# 3'
				start = f.location.start - self.feature_size
				end = f.location.start
				new_features.append(FeatureProperties(start, end, 1, gene_name + "__3prime"))
				
				# 5'
				start = f.location.end
				end = f.location.start + self.feature_size
				new_features.append(FeatureProperties(start, end, 1, gene_name + "__5prime"))
			else:
				# 3'
				start = f.location.end
				end = f.location.start + self.feature_size
				new_features.append(FeatureProperties(start, end, -1, gene_name + "__3prime"))
				
				# 5'
				start = f.location.start - self.feature_size
				end = f.location.start
				new_features.append(FeatureProperties(start, end, -1, gene_name + "__5prime"))
				
				xxxx check the coordinates are within the bounds xxxx
				XXXX check the start and end of the feautres is okay
				XXX complement the feature
				
		return new_features
	
	
	def construct_file(self, filename):
		with open(filename, 'w') as emblfile:
			emblfile.write(self.header())
			
			for f in self.create_3_5_prime_features():
				emblfile.write(self.construct_feature(f))
				
		return self
		
	def feature_to_gene_name(self, feature):
		gene_name_val = 'unknown'
		if "gene" in feature.qualifiers:
			gene_name_val = feature.qualifiers["gene"][0]
		return gene_name_val

	def header(self):
		return """ID   ABC; SV 1; circular; genomic DNA; STD; PRO; {length} BP.
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..{length}
FT                   /organism="Bacteria"
""".format(length=str(self.genome_length))


	def construct_feature(self, feature):
		return """FT   CDS             {window_start}..{window_end}
FT                   /gene="{gene_name}"
FT                   /locus_tag="{gene_name}"
FT                   /product="product"
""".format(gene_name=feature.gene_name(), window_start=str(window.start + 1), window_end=str(window.end))
