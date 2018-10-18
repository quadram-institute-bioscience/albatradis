# This container will install Bio-Tradis from master
#
FROM debian:testing

RUN apt-get update -qq

RUN apt-get install -y sudo smalt samtools tabix make wget unzip zlib1g-dev cpanminus gcc bzip2 libncurses5-dev libncursesw5-dev libssl-dev r-base git python3 python3-setuptools python3-biopython python3-pip 

# Bio-TraDIS toolkit
RUN cpanm -f Bio::Tradis
RUN sudo Rscript -e "source('http://bioconductor.org/biocLite.R')" -e "biocLite(c('edgeR','getopt', 'MASS'))"

# AlbaTraDIS
RUN pip3 install cython
RUN pip3 install git+git://github.com/quadram-institute-bioscience/albatradis.git
