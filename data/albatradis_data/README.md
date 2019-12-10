# Example data for AlbaTraDIS

This folder contains some example data to run the albatradis workflow. The data is from the Triclosan dataset presented in [Yasir et al. 2019](https://www.biorxiv.org/content/10.1101/596833v1). 

## Input files
List of files: 
* reference_BW25113.embl - E-coli BW25113 Genome annotation 
* controlLBrep1.insert_site_plot.gz controlLBrep1.insert_site_plot.gz - LB only control
* 025mgLTricRep1.insert_site_plot.gz 025mgLTricRep2.insert_site_plot.gz - condition plot files (0.25mg/L triclosan concentration)

## Commands

### One condition
This command can be run from this working directory, without any changes, assuming AlbaTraDIS has been successfully installed. It provides a single condition and a single control, which are the minimum required for AlabTraDIS to run. It is a good way to see what the output looks like and to test if the installation has been performed correctly.
```
albatradis -v -a reference_BW25113.embl 025mgLTricRep1.insert_site_plot.gz 025mgLTricRep2.insert_site_plot.gz controlLBrep1.insert_site_plot.gz controlLBrep2.insert_site_plot.gz
```
It takes about 5 minutes to run on a standard laptop, and it is set to output in verbose mode (-v) so that you can see whats happening. 


Alternatively if you dont have AlbaTraDIS installed, you can copy and paste the following docker command (must have docker installed), which will download an installation of AlbaTraDIS and run it without any changes.

```
docker run --rm -it -v $PWD/albatradis_output:/work quadraminstitute/albatradis:latest albatradis -v -a /albatradis/data/albatradis_data/reference_BW25113.embl /albatradis/data/albatradis_data/025mgLTricRep1.insert_site_plot.gz /albatradis/data/albatradis_data/025mgLTricRep2.insert_site_plot.gz /albatradis/data/albatradis_data/controlLBrep1.insert_site_plot.gz /albatradis/data/albatradis_data/controlLBrep2.insert_site_plot.gz
```

### Minimal toy dataset
This is the command to run a minimal toy dataset. It is good for verifying the functionality is working as expected. The results will not be of any significance because there is only a single control and condition.

```
albatradis -v -a reference_BW25113_short.embl 025mgLTricRep1.insert_site_plot_short.gz  controlLBrep1.insert_site_plot_short.gz 
```

Alternatively if you dont have AlbaTraDIS installed, you can copy and paste the following docker command (must have docker installed), which will download an installation of AlbaTraDIS and run it without any changes.
```
docker run --rm -it -v $PWD/albatradis_output:/work quadraminstitute/albatradis:latest albatradis -v -a /albatradis/data/albatradis_data/reference_BW25113_short.embl /albatradis/data/025mgLTricRep1.insert_site_plot_short.gz  /albatradis/data/controlLBrep1.insert_site_plot_short.gz 
```


