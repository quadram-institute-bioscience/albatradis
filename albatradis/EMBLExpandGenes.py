''' Given an annotation file, take each gene, and create a new feature at the start and end to capture promotors'''
from albatradis.EMBLReader import EMBLReader
from albatradis.EMBLSequence import EMBLSequence

class FeatureProperties:
	def __init__(self, 	start, end, direction, gene_name, locus_tag, product):
		self.start = start
		self.end = end
		self.direction = direction
		self.gene_name = gene_name
		self.locus_tag = locus_tag
		self.product = product

class EMBLExpandGenes:
	def __init__(self, embl_file, feature_size):
		self.embl_file = embl_file
		self.feature_size = feature_size
		self.er = EMBLReader(self.embl_file)
		self.features = self.er.features
		self.genome_length = self.er.genome_length
	
	def create_3_5_prime_features(self):
		new_features = []
		for feature in self.features:
			gene_name = self.feature_to_gene_name(feature)
			locus_tag = self.feature_to_locus_tag(feature)
			product = self.feature_to_product(feature)
			
			# The gene itself
			new_features.append(FeatureProperties(feature.location.start, feature.location.end, feature.strand, gene_name, locus_tag, product))
			
			# forward direction
			if feature.strand == 1:
				new_features.append(self.construct_start_feature(feature, gene_name, "__5prime",locus_tag, product))
				new_features.append(self.construct_end_feature(feature, gene_name, "__3prime",locus_tag, product))
			else:
				new_features.append(self.construct_end_feature(feature, gene_name, "__5prime", locus_tag, product))
				new_features.append(self.construct_start_feature(feature, gene_name, "__3prime",locus_tag, product))
					
		return new_features
		
	def construct_end_feature(self, feature, gene_name, suffix, locus_tag, product):
		start = feature.location.end
		end = feature.location.end + self.feature_size
		
		if end > self.genome_length:
			end = self.genome_length
			
		if start >= end or end-start < 10:
			return None
		
		return FeatureProperties(start, end, feature.strand, gene_name + suffix, locus_tag + suffix, product)
		
	def construct_start_feature(self, feature, gene_name, suffix, locus_tag, product):
		start = feature.location.start - self.feature_size
		end = feature.location.start
		
		if start <1: 
			start = 1
		if start >= end or end-start < 10:
			return None
		return FeatureProperties(start, end, feature.strand, gene_name + suffix, locus_tag + suffix, product)
	
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
					
			emblfile.write(EMBLSequence(str(self.er.record.seq)).format())
		return self
		
	def feature_to_gene_name(self, feature):
		gene_name_val = str(feature.location.start) + "_" + str(feature.location.end)
		if "gene" in feature.qualifiers:
			gene_name_val = feature.qualifiers["gene"][0]
		return gene_name_val

	def feature_to_product(self, feature):
		product_val = str(feature.location.start) + "_" + str(feature.location.end)
		if "product" in feature.qualifiers:
			product_val1 = feature.qualifiers["product"][0]
			product_val = product_val1.replace(",", " and ")
		return product_val

	def feature_to_locus_tag(self, feature):
		locus_tag_val = str(feature.location.start) + "_" + str(feature.location.end)
		if "locus_tag" in feature.qualifiers:
			locus_tag_val = feature.qualifiers["locus_tag"][0]
		return locus_tag_val

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
FT                   /locus_tag="{locus_tag}"
FT                   /product="{product}"
""".format(gene_name=feature.gene_name, window_start=str(feature.start +1), window_end=str(feature.end), locus_tag=feature.locus_tag, product=feature.product)

	def construct_feature_reverse(self, feature):
		return """FT   CDS             complement({window_start}..{window_end})
FT                   /gene="{gene_name}"
FT                   /locus_tag="{locus_tag}"
FT                   /product="{product}"
""".format(gene_name=feature.gene_name, window_start=str(feature.start +1), window_end=str(feature.end), locus_tag=feature.locus_tag, product=feature.product)



