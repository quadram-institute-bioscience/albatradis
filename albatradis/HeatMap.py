import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

class HeatMap:
	def __init__(self, reports_to_gene_logfc, gene_names, outputfile):
		self.reports_to_gene_logfc       = reports_to_gene_logfc
		self.gene_names = gene_names
		self.outputfile = outputfile
		
	def clean_filenames(self):
		cleaned_names = []
		
		filenames = self.reports_to_gene_logfc.keys()
		common_prefix = os.path.commonpath(filenames)
		
		for n in self.reports_to_gene_logfc.keys():
			updated_name = n
			if len(common_prefix) > 1:
				updated_name = updated_name.replace(common_prefix,"")
				updated_name = updated_name.strip("/")
			updated_name = updated_name.replace("/gene_report.csv","")
			updated_name = updated_name.replace(".csv","")
			updated_name = updated_name.replace("_"," ")

			cleaned_names.append(updated_name)
			
		return cleaned_names

	def create_heat_map(self):
		int_values = [[int(y) for y in x] for x in self.reports_to_gene_logfc.values()]
		
		dims = (16.54, 23.4)
		fig, ax = plt.subplots(figsize = dims)
		data = pd.DataFrame(int_values, columns = self.gene_names, index=self.clean_filenames())

		heatmap = sns.heatmap(data.T,ax = ax, cmap="RdBu_r", center=0, yticklabels= False)
		
		ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
		plt.tight_layout()
		
		plt.savefig(self.outputfile , dpi=100)
		plt.close()
		return True
