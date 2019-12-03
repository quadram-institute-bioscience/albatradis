# AlbaTraDIS
[![Build Status](https://travis-ci.org/quadram-institute-bioscience/albatradis.svg?branch=master)](https://travis-ci.org/quadram-institute-bioscience/albatradis)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/quadram-institute-bioscience/albatradis/blob/master/LICENSE)
[![Docker Build Status](https://img.shields.io/docker/build/quadraminstitute/albatradis.svg)](https://hub.docker.com/r/quadraminstitute/albatradis)
[![Docker Pulls](https://img.shields.io/docker/pulls/quadraminstitute/albatradis.svg)](https://hub.docker.com/r/quadraminstitute/albatradis)  

## Contents
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Conda](#conda)
    * [Docker](#docker)
	* [Ubuntu/Debian](#ubuntudebian)
  * [Usage](#usage)
  * [License](#license)
  * [Feedback/Issues](#feedbackissues)
  * [Citation](#citation)

## Introduction
AlbaTraDIS is a software application for performing rapid large-scale comparative analysis of TraDIS experiments whilst also predicting the impact of inserts on nearby genes. It allows for experiements with multiple conditions to be easily analysed using statistical methods developed in the Bio-TraDIS toolkit.

## Installation
The software in this repository is straightforward to install, however the Bio-TraDIS toolkit upon which it depends is more complex. If you just want to quickly try out the software please try Docker. If you wish to install it, please use Conda, an finally if you are brave, use the Ubuntu/Debian instructions.

### Conda
[![Anaconda-Server Badge](https://anaconda.org/bioconda/albatradis/badges/latest_release_date.svg)](https://anaconda.org/bioconda/albatradis)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/albatradis/badges/platforms.svg)](https://anaconda.org/bioconda/albatradis)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/albatradis/badges/downloads.svg)](https://anaconda.org/bioconda/albatradis)
```
conda install -c conda-forge -c bioconda albatradis
```

### Docker
Install [Docker](https://www.docker.com/).  There is a docker container which gets automatically built from the latest version of AlbaTraDIS. To install it:

```
docker pull quadraminstitute/albatradis
```

To use it you would use a command such as this (substituting in your filename/directories), using the example file in this respository:
```
docker run --rm -it -v /path/to/example_data:/example_data quadraminstitute/albatradis albatradis xxxxx
```

### Ubuntu/Debian
To install AlbaTraDIS on Ubuntu or Debian run:
```
sudo apt-get update -qq && sudo apt-get install -y bioperl bwa bzip2 cpanminus gcc git libgd-gd2-perl libncurses5-dev libncursesw5-dev libssl-dev libxml-libxml-perl make python3 python3-biopython python3-pip python3-setuptools r-base samtools smalt tabix unzip wget zlib1g-dev
sudo Rscript -e "install.packages('BiocManager')" -e "BiocManager::install()" -e "BiocManager::install(c('edgeR','getopt', 'MASS'))"
cpanm Bio::Tradis

pip3 install cython
pip3 install albatradis
export PATH=${HOME}/.local/bin:$PATH
```
These instructions were tested on Ubuntu 18.04.

## Usage

### albatradis

This is the main script for the application. 


```
usage: albatradis [options] EMBLfile condition_plotfiles control_plotfiles

Tradis analysis

positional arguments:
  emblfile              Annotation file in EMBL format
  plotfiles             Input plot files (optionally gzipped). There must be
                        an equal number of condition and control files

optional arguments:
  -h, --help            show this help message and exit
  --span_gaps SPAN_GAPS, -s SPAN_GAPS
                        Span a gap if it is this multiple of a window size
                        (default: 1)
  --iterations ITERATIONS, -i ITERATIONS
                        No. of times to rescan (default: 1)
  --minimum_block MINIMUM_BLOCK, -b MINIMUM_BLOCK
                        Minimum number of reads which must be in 1 block in
                        comparison (default: 100)
  --minimum_logfc MINIMUM_LOGFC, -f MINIMUM_LOGFC
                        Minimum log fold change +/- (default: 1)
  --minimum_logcpm MINIMUM_LOGCPM, -c MINIMUM_LOGCPM
                        Minimum log counts per million +/- (default: 8.0)
  --minimum_threshold MINIMUM_THRESHOLD, -m MINIMUM_THRESHOLD
                        Only include insert sites with this number or greater
                        insertions (default: 5)
  --minimum_proportion_insertions MINIMUM_PROPORTION_INSERTIONS, -d MINIMUM_PROPORTION_INSERTIONS
                        If the proportion of insertions is too low compared to
                        control, dont call decreased insertions below this
                        level (default: 0.1)
  --dont_normalise_plots, -n
                        Dont normalise input plots (default: False)
  --prefix PREFIX, -o PREFIX
                        Output directory prefix (default: output)
  --pvalue PVALUE, -p PVALUE
                        Dont report anything above this p-value (default:
                        0.05)
  --qvalue QVALUE, -q QVALUE
                        Dont report anything above this q-value (default:
                        0.05)
  --strict_signal, -g   A result must be present in the combined plots to be
                        returned (default: False)
  --use_annotation, -a  Use the reference annotation rather than a sliding
                        window (default: False)
  --prime_feature_size PRIME_FEATURE_SIZE, -z PRIME_FEATURE_SIZE
                        Feature size when adding 5/3 prime block when
                        --use_annotation (default: 198)
  --window_interval WINDOW_INTERVAL, -l WINDOW_INTERVAL
                        Window interval (default: 25)
  --window_size WINDOW_SIZE, -w WINDOW_SIZE
                        Window size (default: 100)
  --verbose, -v         Print out more information about the analysis while it
                        runs (default: False)
  --debug               Turn on debugging (default: False)
  --version             show program's version number and exit
```

#### Positional arguments

__emblfile__: This is an annotated reference genome in EMBL format. It can be downloaded from the EBI website.

__plot_files__: These are insert site plot files, generated using ```bacteria_tradis``` script. The same reference genome must be used in all cases and must match the emblfile. 
Conditions are provided first, followed by controls. The number of conditions must match the number of controls, with a minimum of 1 of each. Ideally you need 2 or more of each.
 
#### Optional arguments

__span_gaps__: When blocks of significat change in insertions are detected they can be fragmented, possibly being at the start and end of a gene and missing in the middle. This option allows you to span these gaps to form more contigous blocks, giving neater results. If you set this too high, then different distinct mechanisms will be merged together, giving you erroneous results.

__iterations__: You can iteratively look for the highest signals, identify them, report them, then mask them out, and start again. This allows you to progressively identify weaker signals which may be overwhelmed ordinarily. There is no automatic stop, so if you do too many iterations, you will increase the number of false positives.

__minimum_block__: This is the minimum number of reads which must be in a block/gene to be considered. If you do a scatter plot experimental variation, low abundance equates to greater variability. This hard minimum threshold eliminates a great deal of noise.

__minimum_logfc__: The minimum log fold change in insertion sites between conditions to consider as significant. It must be an integer. 

__minimum_logcpm__: The minimum log counts per million to consider. It must be an integer and this is approximately equivalent to relative abundance.

__minimum_threshold__: There can be random bits of DNA which makes it through to sequencing, but they are generally uniformly scattered throughout at low frequency. This controls the minimum threshold. Any insert site with less than this number will be set to zero at the start of the experiment.

__minimum_proportion_insertions__: Often experiments produce different numbers of insertions, its just a natural part of the protocol. If the difference between the condition and controls are too extreme, then the statistics start to break down, particularly when calling decreased numbers of insertions, since the absense of data is probably due to the overall number of insertions/reads rather than something real. This is the minimum proportion allowable.

__dont_normalise_plots__: By default plots are normalised to the number of reads in the largest plot file, helping to make the statistics and plots look more uniform. You can turn this off if you wish.

__prefix__: The prefix of the output directory. You'll probably want to change this each time you run the script.

__pvalue__: Dont report anything above this p-value. You may want to reduce the value further, depending on how adventurous you are. If you set it too high you will get more erroneous results, if you set it too low you may not get any results at all.

__strict_signal__: By default if there is a strong signal in the forward or reverse directions above the thresholds, it get reported. You can make this even stricter, by also requiring that the combined data must also have a significant signal. This will reduce the number of genes identified, however may also reduce erroneous signals. 

__use_annotation__: By default the software uses a sliding window to identify significant signals. You can choose instead to use the annotated genes passed in via the input EMBL file. Each gene also has a 5' and 3' feature added, so that you dont miss signals in the intergenic regions and also to identify up and down regulation. 

__prime_feature_size__: When you --use_annotation 5' and 3' features are created around each gene. This controls the size of those features in bases.

__window_interval__: The number of bases to move along when using sliding windows. Ensure it is less than or equal to the windows_size otherwise you will miss parts of the genome. Ideally it should be a maximum of half the window_size.

__window_size__: The size of the sliding window in bases. If you set this too high you will only get very strong signals so will miss quite a bit. If you set this too low you will get a huge amount of false positives due to the natural experimental variation. The window size should be about 10 times the average insertion density, so if there are insertions every 10 bases, the window size should be 100 bases.


#### Output files

__annotation.embl__: This file contains a modified version of the input annotation file in EMBL format. For example, it adds 5' and 3' features to each gene. It can be visualised in Artemis.

__gene_report.csv__: This is the main results file in the output of the script. It contains a tab delimited spreadsheet detailing genes identified as being interesting, in the below format. The first column is the gene name derived from the input annotation file. If a signal is identified in an intergenic region, the start and end coordinates are given. The next column is a categorisation of the mechanism, such as up or down regulation, knockout, unclassified etc.... The 3rd and 4th columns are the coordinates of the start and end of the gene, relative to the input annotation file. The MaxLogFC is the maximum log fold change in the signal observed in the gene (or in the 5'/3'). It is rounded to the nearest integer.  If the signal is strongest in a single direction, the max log fold change in that direction is reported rather than the value for the combined analysis. The expression column indicates if the gene is experiencing an increase or decrease in insertions. The direction column indicates which direction the significant insertions where primarily detected in (or no direction if both apply). Finally the last column gives the upstream gene, which is often implicated in the mechanism.

|  Gene | Category | Start | End | MaxLogFC | Expression | Direction | Upstream |
|  --- | --- | --- | --- | --- | --- | --- | --- |
|  zabC | downregulated | 100 | 500 | 1 | increased_insertions | forward | abc |
|  yxxY | upregulated | 135 | 234 | 1 | decreased_insertions | reverse | efg |

__regulated_gene_report.csv__: This is a filtered version of the gene_report.csv file but only includes genes which are identified as upregulated or downregulated. If no genes are identified with this pattern, the file will not be created.  This is really only useful if the experiements included a promotor.

__combined.csv__: This comma delimited spreadsheet is the output of the Bio-TraDIS toolkit, with additional essentiality categoristations, and an example is listed below. This is the raw data from which the gene_report.csv is derived. It lists each gene or sliding window, and optionally the corresponding 5' and 3' features for a gene. The first 2 columns list the names of the gene or give the coordinates of the sliding window. The 3rd column lists the annotated function of the gene (if available in the annotation file). The numerical columns are derived from EdgeR. The 4th column gives the log fold change between the conditions and the controls. The 5th column gives the log counts per million, which can be thought of as relative abundance. The final column indicate how the essentiality has changed between the conditions and the controls, so a gene can always be non-essential in both the controls and the conditions or essential in all cases. More interestingly though is where there is a change in essentiality between the control and the conditions, indicating a large mechanistic change. 

| locus_tag | gene_name | function | logFC | logCPM | PValue | q.value | Essentiality |
| --- | --- | --- | --- | --- | --- | --- | --- |
| thrL | thrL | product | -0.4327 | 4.1269 | 0.5477 | 0.8177 | always_nonessential |
| thrL__5prime | thrL__5prime | product | -0.1208 | 4.5885 | 0.8555 | 0.9521 | always_nonessential |
| thrL__3prime | thrL__3prime | product | 1.0268 | 4.9723 | 0.1227 | 0.4258 | always_nonessential |

__forward.csv__: This is identical to the combined.csv file, except only insertions in the forward direction were considered during the analysis.

__reverse.csv__: This is identical to the combined.csv file, except only insertions in the reverse direction were considered during the analysis.

__combined.plot__: This is the log fold change of each gene or sliding window, in a User plot format suitable for viewing in Artemis. It consists of 2 space delimited integers on each line, where a line corresponds to a base in the reference genome. A positive integer means there has been an increase in insertions, and a negative integer means there has been a decrease in insertions.

__forward.plot__: This is identical to the combined.plot file, except only insertions in the forward direction were considered during the analysis. 

__reverse.plot__: This is identical to the combined.plot file, except only insertions in the reverse direction were considered during the analysis. 

### albatradis-presence_absence

After you have run the albatradis script, it produces gene_report files. This script performs comparitive analysis and outputs heatmaps, combined spreadsheets, figures and trees (dendrograms).

```
usage: albatradis-presence_absence [options] EMBLfile gene_reports

Take in gene report files and produce a heatmap

positional arguments:
  emblfile              Annotation file in EMBL format
  genereports           Gene report spreadsheets

optional arguments:
  -h, --help            show this help message and exit
  --prefix PREFIX, -o PREFIX
                        Output directory prefix (default: output)
  --verbose, -v         Print out more information about the analysis while it
                        runs (default: False)
  --debug               Turn on debugging (default: False)
  --version             show program's version number and exit
```

#### Positional arguments

__emblfile__: This is an annotated reference genome in EMBL format. It can be downloaded from the EBI website.

__genereports__: One of the outputs of the albatradis script is a gene_report.csv file. You will have 1 of these for each condition, and so providing all of them here will allow for the comparison of conditions. Its probably best to add the name of your condition into the file name to make it easier to identify in the output.

#### Optional arguments

__prefix__: This is the output directory prefix and there are a number of output files.


#### Output files
__all_logfc.csv__: A tab delimited spreadsheet containing the log fold change integer value for every gene against every condition. It is usually a huge big sparse matrix, so its best to process it with a script.

| Sample | aaaA | bbbB | cccC | dddD |
| ------ | ---- | ---- | ---- | ---- |
| Cond1  | 1    | -1   | 0    | 0    |
| Cond2  | 0    | 1    | 16   | 0    |
| Cond3  | 0    | 8    | 0    | 0    |
| Cond4  | 1    | 0    | 9    | 0    |

__filtered_logfc.csv__: Identical to the all_logfc.csv file, except genes only genes with a significant signal in at least 1 condition are kept. This is much easier to look at manually.

__full_heatmap.png__: A visual representation of the data in all_logfc.csv in PNG format.

__filtered_heatmap.png__: A visual representation of the data in filtered_logfc.csv in PNG format.

__distance_matrix_dendrogram.png__: A simple dendrogram  (tree) figure in PNG format of the relationships between the conditions, based on the number of shared genes with significant signals. 

__nj_newick_tree.tre__: A neighbour joining tree in Newick format, created from a distance matrix, based on the number of shared genes with significant signals. This tree can be visualised in FigTree or any number of different applications.

__logfc.dot__: A graph representation of the shared genes between conditions in DOT format. This is a standard graphing format, with substantial support API supoort. It can be interactively visualised with Gephi. There is 1 node for each condition and signficiatn gene. The edges represent where a gene is found in a condition, linking the two nodes. This then nicely shows the network of shared mechanisms of action.

__union_gene_report.csv__: A comma separated spreadsheet in the same format as the gene_report.csv file, consisting of a union of all of the input files. A gene is represented by 1 row.

__intersection_gene_report.csv__: A comma separated spreadsheet in the same format as the gene_report.csv file, consisting of the intersection of all of the input files. So only genes which are found in every condition (common modes of action) are in the file. A gene is represented by 1 row.

### albatradis-gene_reports
This is a helper script that you may never need as the functionality is used within the albatradis-presence_absence script. You can take multiple gene_report.csv files and perform set operations on them. It is useful if you know that a few conditions should be merged together as the mechanisms are identical.

```
usage: albatradis-gene_reports [options] gene_report1.csv gene_report2.csv ...

Manipulate gene_report.csv files, such as performing set operations

positional arguments:
  genereports           Gene report spreadsheets

optional arguments:
  -h, --help            show this help message and exit
  --prefix PREFIX, -o PREFIX
                        Output directory prefix (default: output)
  --verbose, -v         Print out more information about the analysis while it
                        runs (default: False)
  --debug               Turn on debugging (default: False)
  --version             show program's version number and exit
```

#### Positional arguments

__genereports__: One of the outputs of the albatradis script is a gene_report.csv file. You will have 1 of these for each condition.

#### Optional arguments

__prefix__: This is the output directory prefix and there are a number of output files.

#### Output files

__union_gene_report.csv__: A comma separated spreadsheet in the same format as the gene_report.csv file, consisting of a union of all of the input files. A gene is represented by 1 row.

__intersection_gene_report.csv__: A comma separated spreadsheet in the same format as the gene_report.csv file, consisting of the intersection of all of the input files. So only genes which are found in every condition (common modes of action) are in the file. A gene is represented by 1 row.


### albatradis-scatterplot

This script produces scatterplots of your input data plotted against itself and the controls. It is useful as a QC metric to see if the data is biased. Basically you take sliding windows, count the number of reads in each window, then plot those values against the other condition and against the others. This is on a log scale and the outliers are the interesting points.

```
usage: albatradis-scatterplot [options] --control control1.plot --control control2.plot --condition condition1.plot --condition condition2.plot

Create scatter plot of controls vs conditions

optional arguments:
  -h, --help            show this help message and exit
  --control CONTROL, -c CONTROL
                        control files (use 2 or more) (default: None)
  --condition CONDITION, -d CONDITION
                        condition files (use 2 or more) (default: None)
  --window_size WINDOW_SIZE, -w WINDOW_SIZE
                        Window size (default: 50)
  --outputfile OUTPUTFILE, -o OUTPUTFILE
                        Output filename prefix (default: scatter)
  --normalise, -n       normalise the files (default: False)
  --verbose, -v         Print out more information while it runs (default:
                        False)
  --debug               Turn on debugging (default: False)
  --version             show program's version number and exit
```

#### Positional arguments

__control___: The insert site plots of the controls, where you must have 2 or more.

__condition__: The insert site plots of the conditions, where you must have 2 or more.

#### Optional arguments

__window_size___: The size of the window in bases. The interval is set to the window_size.

__outputfile__: The output file prefix.

__normalise__: Normalise the input files reads to the input file with the largest number of reads.


#### Output files

__scatter*.png__: Images of the scatterplots on a log scale in PNG format.


### albatradis-annotation
Take in an EMBL file and add flanking 3 prime and 5 prime annotation. It is used as part of the albatradis --use_annotation feature, so you may not need it, as the annotated file is saved in the output directory.  

```
usage: albatradis-annotation [options] EMBLfile

Take in an EMBL file and add flanking 3 prime and 5 prime annotation

positional arguments:
  emblfile              Annotation file in EMBL format

optional arguments:
  -h, --help            show this help message and exit
  --feature_size FEATURE_SIZE, -s FEATURE_SIZE
                        Feature size (default: 198)
  --outputfile OUTPUTFILE, -o OUTPUTFILE
                        Output file (default: output.embl)
  --verbose, -v         Print out more information about the analysis while it
                        runs (default: False)
  --debug               Turn on debugging (default: False)
  --version             show program's version number and exit
```

#### Positional arguments

__emblfile__: This is an annotated reference genome in EMBL format. It can be downloaded from the EBI website.

#### Optional arguments

__feature_size__: When you --use_annotation 5' and 3' features are created around each gene. This controls the size of those features in bases.

__outputfile__: The name of the output file

#### Output files

__output.embl__: The original EMBL file, plus annotated 5' and 3' features, to give another EMBL file, including the reference genome sequence.


### albatradis-artemis_project
Sometimes you want to view the insert site plots in Artemis. It can be quite a manual task to open up different replicates and combinations. This script will generate a project.properties file from a spreadsheet which gets automatically loaded by Artemis (from the current working directory). This then makes it quicker to view multiple different insert site plots.

```
usage: albatradis-artemis_project [options] reference experiments_metadata.csv

Create an artemis project file

positional arguments:
  reference             reference EMBL file
  experiments_metadata  experiments metadata spreadsheet

optional arguments:
  -h, --help            show this help message and exit
  --control CONTROL, -c CONTROL
                        control files (can use multiple times) (default: None)
  --outputfile OUTPUTFILE, -o OUTPUTFILE
                        Output filename (default: project.properties)
  --verbose, -v         Print out more information while it runs (default:
                        False)
  --debug               Turn on debugging (default: False)
  --version             show program's version number and exit
```
#### Positional arguments

__reference__: This is an annotated reference genome in EMBL format. It can be downloaded from the EBI website.

__experiments_metadata__: A comma separated spreadsheet in the format of "Drug,Pathway,DetailedPathway,Impact,MIC,Induction,Rep1,Rep2".

#### Optional arguments

__control__: Path to the control files.

__outputfile__: The name of the Artemis project file. If you change this, then Artemis wont work.

#### Output files

__project.properties__: The Artemis project file.

## License
AlbaTraDIS is free software, licensed under [GPLv3](https://raw.githubusercontent.com/quadram-institute-bioscience/albatradis/master/VERSION/LICENSE).

## Feedback/Issues
Please report any issues or to provide feedback please go to the [issues page](https://github.com/quadram-institute-bioscience/albatradis/issues). If you make improvements to the software, please send us the changes though a [pull request](https://github.com/quadram-institute-bioscience/albatradis/pulls) so that the whole community may benefit from your work.

## Citation
["AlbaTraDIS: Comparative analysis of large datasets from parallel transposon mutagenesis experiments"](https://doi.org/10.1101/593624), Andrew J. Page, Sarah Bastkowski, Muhammad Yasir, A. Keith Turner, Thanh Le Viet, George M. Savva, Mark A. Webber, Ian G. Charles
bioRxiv 593624; doi: https://doi.org/10.1101/593624


