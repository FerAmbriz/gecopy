FROM debian:stable
#-------------------------- Dependencies -------------------------#
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-dev python3-pip r-base git gfortran curl libcurl4-openssl-dev libbz2-dev pkg-config libcairo2-dev libjpeg-dev libgif-dev \
  make gcc g++ zlib1g-dev liblapack-dev libblas-dev libpng-dev liblzma-dev && ln -s /usr/bin/python3 /usr/bin/python && rm -rf /var/lib/apt/lists/* && pip install --upgrade pip setuptools wheel && pip install --no-cache-dir CNVkit

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('panelcn.mops', 'Rhtslib', 'Rsamtools', 'GenomicAlignments', 'DNAcopy'), ask=FALSE)"
RUN R -e "install.packages('ExomeDepth')"
COPY scr /usr/bin/
