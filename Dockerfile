FROM continuumio/miniconda
ARG app_version

MAINTAINER Andrew McPherson <andrew.mcpherson@gmail.com>

USER root
RUN conda config --add channels https://conda.anaconda.org/dranew && conda config --add channels https://conda.anaconda.org/shahcompbio && conda config --add channels bioconda
RUN conda install defuse==$app_version
RUN conda install openssl=1.0
RUN mkdir -p /root/.config/matplotlib
RUN echo "backend : Agg" > /root/.config/matplotlib/matplotlibrc

ENV NAME defuse

ENTRYPOINT ["bash"]

