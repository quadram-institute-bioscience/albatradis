class EMBLGenerator:
	def __init__(self, windows, genome_length):
		self.windows = windows
		self.genome_length = genome_length
		
	def construct_file(self, filename):
		with open(filename, 'w') as emblfile:
			emblfile.write(self.header())
			
			for w in self.windows:
				emblfile.write(self.window_feature(w))
				
		return self
		

	def header(self):
		return """ID   ABC; SV 1; circular; genomic DNA; STD; PRO; {length} BP.
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..{length}
FT                   /organism="Bacteria"
""".format(length=str(self.genome_length))


	def window_feature(self, window):
		return """FT   CDS             {window_start}..{window_end}
FT                   /gene="{gene_name}"
FT                   /locus_tag="{gene_name}"
FT                   /product="product"
""".format(gene_name=window.name_start_one(), window_start=str(window.start + 1), window_end=str(window.end))

