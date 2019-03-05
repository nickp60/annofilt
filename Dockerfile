FROM continuumio/miniconda3
MAINTAINER Nick Waters <nickp60@gmail.com>
RUN conda install -c BioBuilds blast biopython
RUN git clone https://github.com/nickp60/annofilt.git ######
# ADD . annofilt
WORKDIR annofilt
RUN python setup.py install
ENTRYPOINT [ "annofilt" ]
