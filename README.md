# scIsoSeq
This is a [nextflow](https://github.com/nextflow-io/nextflow) pipeline for analyzing multiplexed single-cell IsoSeq (PacBio long-read RNA sequencing). This is configured particularly for data from Smart-Seq2 cDNA libraries with additional barcodes added for multiplexing on a PacBio system. Currently, it assumes that the input is indexed CCS reads.

## Setup
### Requirements
This pipeline requires [anaconda](https://anaconda.org/). The required dependencies will then be installed by nextflow into a conda virtual environment

### Configuration
The current pipeline is designed to run with a SLURM scheduler. 

To make it run for your group edit the `projectName` & `fasta` parameters in the the mendel.config file to fit to your group project and needs. To make it run on another infrastructure simply add a new nextflow config file in the conf folder and source via the nextflow.config file. See [here](https://www.nextflow.io/docs/latest/config.html) for more information. 

## Workflow of the pipeline
2. demuxing and primer removal (lima)
3. refine reads (isoseq3)
4. merge samples (optional)
3. cluster reads 
4. polish reads
5. align transcripts (minimap2)
6. collapse redudant transcripts (Cupcake ToFU)
7. Compare with annotation (SQANTI2)
8. Filter artifacts (SQANTI2) (optionally with short-read junction data)



## Run the pipeline
To run the pipeline simply run e.g:

```bash
nextflow run_scIsoSeq.nf --input "/path/to/input" --output "/path/to/output"
