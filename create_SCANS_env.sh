#!/bin/bash

# Environment to create
ENV_PATH="./SCANS_ENV"

# create environment with principal packages
conda create -y -p $ENV_PATH \
 bioconda::feelnc=0.2 \
 bioconda::liftoff=1.6.3 \
 bioconda::bedtools=2.30.0 \
 bioconda::ucsc-genepredtobed=377 \
 bioconda::ucsc-gtftogenepred=377 \
 perl \
 bioconda::perl-bio-featureio

# activation
conda activate $ENV_PATH

# install R packages
conda install -y \
conda-forge::r-base=4.1.3 \
  bioconda::bioconductor-rtracklayer=1.54.0

conda install -y \
  conda-forge::r-tidyverse=1.3.2 \
  conda-forge::r-pbapply=1.7_0 \
  conda-forge::r-stringi=1.7.12 \
  conda-forge::r-stringr=1.5.0 \
  conda-forge::r-tidyr=1.3.0 \
  conda-forge::r-ggplot2=3.4.2 \
  conda-forge::r-ggpubr=0.6.0 \
  conda-forge::r-upsetr=1.4.0 \

echo "Conda environment installed: $ENV_PATH"
echo "To activate SCANS environment: conda activate $ENV_PATH"
