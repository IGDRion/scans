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
  bioconda::perl-bio-featureio \
  conda-forge::r-base \
  bioconda::bioconductor-rtracklayer \
  conda-forge::r-tidyverse \
  conda-forge::r-pbapply \
  conda-forge::r-stringi \
  conda-forge::r-stringr \
  conda-forge::r-tidyr \
  conda-forge::r-ggplot2 \
  conda-forge::r-ggpubr \
  conda-forge::r-upsetr \
&& echo "Conda environment installed: $ENV_PATH" \
&& echo "To activate SCANS environment: conda activate $ENV_PATH"
