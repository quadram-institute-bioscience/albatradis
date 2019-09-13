This folder contains some example data to run the albatradis workflow. The data is from the Triclosan dataset presented in Yasir et al. 2019 (https://www.biorxiv.org/content/10.1101/596833v1). 
List of files: 
	reference_BW25113.embl - E-coli BW25113 Genome annotation 
	controlLBrep1.insert_site_plot.gz controlLBrep1.insert_site_plot.gz - LB only control
	025mgLTricRep1.insert_site_plot.gz 025mgLTricRep2.insert_site_plot.gz - condition plot files (0.25mg/L triclosan concentration)
Example command to run albatradis:

albatradis -v -a reference_BW25113.embl 025mgLTricRep1.insert_site_plot.gz 025mgLTricRep2.insert_site_plot.gz controlLBrep1.insert_site_plot.gz controlLBrep1.insert_site_plot.gz

This runs albatradis script using the annotation and verbose mode using example data. Please run albatradis -h for more options.

