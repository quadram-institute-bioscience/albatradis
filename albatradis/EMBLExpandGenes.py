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
		er = EMBLReader(self.embl_file)
		self.features = er.read_annotation_features()
		self.genome_length = er.genome_length
	
	def create_3_5_prime_features(self):
		new_features = []
		for feature in self.features:
			gene_name = self.feature_to_gene_name(feature)
			
			# The gene itself
			new_features.append(FeatureProperties(feature.location.start, feature.location.end, feature.strand, gene_name))
			
			# forward direction
			if feature.strand == 1:
				new_features.append(self.construct_start_feature(feature, gene_name, "__3prime"))
				new_features.append(self.construct_end_feature(feature, gene_name, "__5prime"))
			else:
				new_features.append(self.construct_end_feature(feature, gene_name, "__3prime"))
				new_features.append(self.construct_start_feature(feature, gene_name, "__5prime"))
					
		return new_features
		
	def construct_end_feature(self,feature, gene_name, suffix):
		start = feature.location.end
		end = feature.location.end + self.feature_size
		
		if end > self.genome_length:
			end = self.genome_length
			
		if start >= end or end-start < 10:
			return None
		
		return FeatureProperties(start, end, feature.strand, gene_name + suffix)
		
	def construct_start_feature(self,feature, gene_name, suffix):
		start = feature.location.start - self.feature_size
		end = feature.location.start
		
		if start <1: 
			start = 1
		if start >= end or end-start < 10:
			return None
		return FeatureProperties(start, end, feature.strand, gene_name + suffix)
	
	def construct_file(self, filename):
		with open(filename, 'w') as emblfile:
			emblfile.write(self.header())
			
			for f in self.create_3_5_prime_features():
				if f == None:
					continue
				if f.direction == 1:
					emblfile.write(self.construct_feature_forward(f))
				else:
					emblfile.write(self.construct_feature_reverse(f))
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


	def construct_feature_forward(self, feature):
		return """FT   CDS             {window_start}..{window_end}
FT                   /gene="{gene_name}"
FT                   /locus_tag="{gene_name}"
FT                   /product="product"
""".format(gene_name=feature.gene_name, window_start=str(feature.start +1), window_end=str(feature.end))

	def construct_feature_reverse(self, feature):
		return """FT   CDS             complement({window_start}..{window_end})
FT                   /gene="{gene_name}"
FT                   /locus_tag="{gene_name}"
FT                   /product="product"
""".format(gene_name=feature.gene_name, window_start=str(feature.start +1), window_end=str(feature.end))

