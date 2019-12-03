FROM sangerpathogens/bio-tradis:latest

RUN apt-get update -qq && apt-get install -y sudo python3 python3-setuptools python3-biopython python3-pip gcc

RUN pip3 install cython
ADD . /albatradis
RUN pip3 install /albatradis
WORKDIR /work
