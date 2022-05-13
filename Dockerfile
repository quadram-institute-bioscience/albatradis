FROM sbastkowski/quatradis:0.8.2

ADD . /albatradis
WORKDIR /albatradis
RUN pip3 install -r /albatradis/requirements.txt
RUN pip3 install .[dev]
RUN pip3 install .

WORKDIR /work
