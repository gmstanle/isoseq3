LABEL author="Geoff Stanley" \
      description="Docker image containing all requirements for scIsoSeq3"

# install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to path (instead of doing 'conda activate')
