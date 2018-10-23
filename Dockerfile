# This container will install Bio-Tradis from master
#
FROM debian:testing
RUN apt-get update -qq
RUN apt-get install -y sudo bio-tradis git python3 python3-setuptools python3-biopython python3-pip 

# AlbaTraDIS
RUN pip3 install cython
RUN pip3 install git+git://github.com/quadram-institute-bioscience/albatradis.git
