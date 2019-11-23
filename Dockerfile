FROM nfcore/base:1.7
LABEL author="Geoff Stanley" \
      description="Docker image containing all requirements for scIsoSeq"


# install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-scIsoSeq/bin:$PATH


