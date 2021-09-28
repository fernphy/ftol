FROM rocker/rstudio:4.1.1

ARG DEBIAN_FRONTEND=noninteractive

############################
### Install APT packages ###
############################

# gawk for taxon-tools
# gcc through libtool for treePL
# cmake, libeigen3-dev for IQTREE
# zlib1g-dev for R package XVector
# libxml2-dev for R package XML
# libudunits2-dev for R package units
# libgdal-dev for R package sf
# libzmq3-dev for R package rzmq -> clustermq
# libmagick++-dev for R package magick -> phytools

# curl, python-dev-is-python3 for hybpiper
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    libnlopt-dev \
    libnlopt0 \
    libcolpack-dev \
    make \
    libomp-dev \
    build-essential \
    autoconf \
    autotools-dev \
    automake \
    libtool \
    cmake \
    libeigen3-dev \
    zlib1g-dev \
    libxml2-dev \
    libudunits2-dev \
    libgdal-dev \
    libzmq3-dev \
    libmagick++-dev \
    mafft \
    ncbi-blast+ \
    fastp \
    time \
    parallel \
    python-dev-is-python3 \
    curl \
    fasttree \
    gawk \
  && apt-get clean

########################
### python libraries ###
########################

# biopython for Hybpiper

RUN curl https://bootstrap.pypa.io/get-pip.py | python \
  && pip install biopython

#############################
### Other custom software ###
#############################

ENV APPS_HOME=/apps
RUN mkdir $APPS_HOME
WORKDIR $APPS_HOME

### treePL ###
RUN git clone https://github.com/blackrim/treePL.git \
  && cd $APPS_HOME/treePL/deps/ \
  && tar xvzf adol-c_git_saved.tar.gz \
  && cd $APPS_HOME/treePL/deps/adol-c/ \
  && ./update_versions.sh \
  && ./configure --with-openmp-flag=-fopenmp --prefix=/usr \
  && make \
  && make install \
  && cd $APPS_HOME/treePL/src \
  && ./configure \
  && make \
  && echo '/usr/lib64' > /etc/ld.so.conf.d/lib64.conf \
  && ldconfig \
  && cp treePL /usr/local/bin

### gnparser ###
WORKDIR $APPS_HOME
ENV APP_NAME=gnparser
ENV VERSION=1.4.0
ENV DEST=$APPS_HOME/$APP_NAME/$VERSION
RUN wget https://github.com/gnames/gnparser/releases/download/v$VERSION/gnparser-v$VERSION-linux.tar.gz \
  && tar xf $APP_NAME-v$VERSION-linux.tar.gz \
  && rm $APP_NAME-v$VERSION-linux.tar.gz \
  && mv "$APP_NAME" /usr/local/bin/

### IQ Tree ###
WORKDIR $APPS_HOME
ENV APP_NAME=IQ-TREE
RUN git clone https://github.com/Cibiv/$APP_NAME.git && \
	cd $APP_NAME && \
	mkdir build && \
	cd build && \
	cmake -DIQTREE_FLAGS=omp .. && \
	make && \
	cp iqtree /usr/local/bin

### trimAL ###
WORKDIR $APPS_HOME
ENV APP_NAME=trimal
RUN git clone https://github.com/scapella/$APP_NAME.git && \
	cd $APP_NAME/source && \
	make && \
	cp trimal /usr/local/bin

### HybPiper ###
# Use my fork until PR gets submitted and merged (or chmod won't work).
WORKDIR $APPS_HOME
ENV APP_NAME=HybPiper
ENV DEST=$APPS_HOME/$APP_NAME
RUN git clone https://github.com/joelnitta/$APP_NAME.git \
  && chmod +x $DEST/*.py
ENV PATH="$PATH:$DEST/"

# Silence parallel citation warning
RUN mkdir /root/.parallel \
  && touch /root/.parallel/will-cite
  
### taxon-tools ###
WORKDIR $APPS_HOME
ENV APP_NAME=taxon-tools
RUN git clone https://github.com/camwebb/$APP_NAME.git && \
	cd $APP_NAME && \
  git checkout 8f8b5e2611b6fdef1998b7878e93e60a9bc7c130 && \
	make check && \
	make install

######################
### conda packages ###
######################

# install miniconda
ENV MINICONDA_VERSION py37_4.8.3
ENV CONDA_DIR $HOME/miniconda3

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-$MINICONDA_VERSION-Linux-x86_64.sh -O ~/miniconda.sh && \
    chmod +x ~/miniconda.sh && \
    ~/miniconda.sh -b -p $CONDA_DIR && \
    rm ~/miniconda.sh

# make non-activate conda commands available
ENV PATH=$CONDA_DIR/bin:$PATH

# make conda activate command available from /bin/bash --login shells
RUN echo ". $CONDA_DIR/etc/profile.d/conda.sh" >> ~/.profile

# make conda activate command available from /bin/bash --interative shells
RUN conda init bash

# build conda environment: kindel

ENV APPNAME kindel

COPY conda/$APPNAME-environment.yml /tmp/

ENV ENV_PREFIX /env/$APPNAME
RUN conda update --name base --channel defaults conda && \
    conda env create --prefix $ENV_PREFIX --file /tmp/$APPNAME-environment.yml --force && \
    conda clean --all --yes

# make shell script to run conda app
# e.g., kindel can be run with `kindel`

RUN echo '#!/bin/bash' >> /usr/local/bin/$APPNAME && \
  echo "source /miniconda3/etc/profile.d/conda.sh" >> /usr/local/bin/$APPNAME && \
  echo "conda activate /env/$APPNAME" >> /usr/local/bin/$APPNAME  && \
  echo "$APPNAME \"\$@\"" >> /usr/local/bin/$APPNAME  && \
  chmod 755 /usr/local/bin/$APPNAME

# build conda environment: bbmap

ENV APPNAME bbmap

COPY conda/$APPNAME-environment.yml /tmp/

ENV ENV_PREFIX /env/$APPNAME
RUN conda update --name base --channel defaults conda && \
    conda env create --prefix $ENV_PREFIX --file /tmp/$APPNAME-environment.yml --force && \
    conda clean --all --yes

# make shell script to run conda app
# NOTE: the actual command within conda to run bbmap is bbmap.sh

RUN echo '#!/bin/bash' >> /usr/local/bin/$APPNAME && \
  echo "source /miniconda3/etc/profile.d/conda.sh" >> /usr/local/bin/$APPNAME && \
  echo "conda activate /env/$APPNAME" >> /usr/local/bin/$APPNAME  && \
  echo "$APPNAME.sh \"\$@\"" >> /usr/local/bin/$APPNAME  && \
  chmod 755 /usr/local/bin/$APPNAME

####################################
### Install R packages with renv ###
####################################

# Create directory for renv project library
RUN mkdir renv

# Modify Rprofile.site so renv uses /renv for project library
RUN echo 'Sys.setenv(RENV_PATHS_LIBRARY = "/renv")' >> /usr/local/lib/R/etc/Rprofile.site

# Initialize a 'dummy' project and restore the renv library.
# Since the library path is specified as above, the library will be restored to /renv
RUN mkdir /tmp/project

COPY ./renv.lock /tmp/project

WORKDIR /tmp/project

# Don't use cache (the symlinks won't work from Rstudio server)
RUN Rscript -e 'install.packages("renv"); renv::consent(provided = TRUE); renv::settings$use.cache(FALSE); renv::init(bare = TRUE); renv::restore()'

WORKDIR /home/rstudio/
