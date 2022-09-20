
################## BASE IMAGE ######################
FROM rocker/tidyverse
################## METADATA ######################
LABEL base_image="rocker/tidyverse"
LABEL version="0"
LABEL software="nf-compare-miRNome"
LABEL software.version="0.2"
LABEL about.summary="Nextflow pipeline to extract and create miRNA and 3'UTR consensus sequences for analysis \n with targetscan and miRmap"
LABEL about.home="https://github.com/Ed-G655/nf-miRNA-SNPs-classify.nf.git"
LABEL about.documentation="https://github.com/Ed-G655/nf-miRNA-SNPs-classify.nf.git"
LABEL about.tags="Genomics"

################## MAINTAINER ######################
MAINTAINER Jose Eduardo Garcia Lopez <jeduardogl655@gmail.com>

################## INSTALLATION ######################

# install pipeline requirements
RUN apt-get update && apt-get install -y --no-install-recommends \
  gcc-multilib \
  apt-utils



# install python3 and pip
RUN apt-get install -y --no-install-recommends \
      python3.5 \
      python3-pip \
      && apt-get clean && \
      rm -rf /var/lib/apt/lists/*

# install python dependecies
RUN pip install numpy
RUN pip install biopython


# install R packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  && install2.r --error \
  pacman
