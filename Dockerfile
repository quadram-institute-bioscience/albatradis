FROM debian:testing
RUN apt-get update -qq && apt-get install -y sudo bio-tradis git python3 python3-setuptools python3-biopython python3-pip cpanminus libncursesw5-dev libssl-dev
RUN cpanm -f IPC::System::Simple DateTime::Locale DateTime Bio::Tradis
RUN Rscript -e "source('http://bioconductor.org/biocLite.R')" -e "biocLite(c('edgeR','getopt', 'MASS'))"

# AlbaTraDIS
RUN pip3 install cython
RUN pip3 install git+git://github.com/quadram-institute-bioscience/albatradis.git
