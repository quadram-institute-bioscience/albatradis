import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

class HeatMap:
	def __init__(self, reports_to_gene_logfc, gene_names, outputfile):
		self.reports_to_gene_logfc       = reports_to_gene_logfc
		self.gene_names = gene_names
		self.outputfile = outputfile
		
	def clean_filenames(self):
		cleaned_names = []
		for n in self.reports_to_gene_logfc.keys():
			cleaned_names.append(n.replace("/gene_report.csv",""))
		return cleaned_names

	def create_heat_map(self):
		int_values = [[int(y) for y in x] for x in self.reports_to_gene_logfc.values()]
		
		a4_dims = (8.27, 11.7)
		fig, ax = plt.subplots(figsize=a4_dims)
		data = pd.DataFrame(int_values, columns = self.gene_names, index=self.clean_filenames())

		heatmap = sns.heatmap(data.T,ax = ax, cmap="RdBu_r", center=0, yticklabels= False)
		plt.savefig(self.outputfile , dpi=100)
		plt.close()
		