FROM nfcore/base:1.7
LABEL author="Geoff Stanley" \
      description="Docker image containing all requirements for scIsoSeq"


# install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-scIsoSeq/bin:$PATH


