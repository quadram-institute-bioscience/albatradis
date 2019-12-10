# Example data for AlbaTraDIS albatradis-presence_absence script

This folder contains some example data to run the albatradis-presence_absence workflow. The data is from the Triclosan dataset presented in Yasir et al. 2019 (https://www.biorxiv.org/content/10.1101/596833v1). 

## Input files
List of files: 
* reference_BW25113.embl - E-coli BW25113 Genome annotation 
* gene_report_0008mgL.csv - gene report for 0.008mg Triclosan concentration
* gene_report_0015mgL.csv - gene report for 0.015mg Triclosan concentration
* gene_report_003mgL.csv - gene report for 0.03mg Triclosan concentration
* gene_report_006mgL.csv - gene report for 0.06mg Triclosan concentration
* gene_report_0125mgL.csv - gene report for 0.125mg Triclosan concentration
* gene_report_025mgL.csv - gene report for 0.25mg Triclosan concentration
* gene_report_05mgL.csv - gene report for 0.5mg Triclosan concentration
* gene_report_1mgL.csv - gene report for 1mg Triclosan concentration


## Commands

This command can be run from this working directory, without any changes, assuming AlbaTraDIS has been successfully installed. It is a good way to see what the output looks like and to test if the installation has been performed correctly.
```
albatradis-presence_absence reference_BW25113.embl gene_report_0008mgL.csv gene_report_0015mgL.csv gene_report_003mgL.csv gene_report_006mgL.csv gene_report_0125mgL.csv gene_report_025mgL.csv gene_report_05mgL.csv gene_report_1mgL.csv
```
This command only takes a few seconds to run on a standard laptop.

Alternatively if you don't have AlbaTraDIS installed, you can copy and paste the following docker command (must have docker installed), which will download an installation of AlbaTraDIS and run it without any changes.
```
docker run --rm -it -v $PWD/data/presence_absence_data:/work quadraminstitute/albatradis:latest albatradis-presence_absence reference_BW25113.embl gene_report_0008mgL.csv gene_report_0015mgL.csv gene_report_003mgL.csv gene_report_006mgL.csv gene_report_0125mgL.csv gene_report_025mgL.csv gene_report_05mgL.csv gene_report_1mgL.csv
```