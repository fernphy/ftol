FROM rocker/r-ver:4.3.1

ARG DEBIAN_FRONTEND=noninteractive

############################
### Install APT packages ###
############################

# gawk for taxon-tools
# gcc through libtool for treePL
# zlib1g-dev for R package XVector
# libxml2-dev for R package XML
# libudunits2-dev for R package units
# libgdal-dev for R package sf
# libzmq3-dev for R package rzmq -> clustermq
# libmagick++-dev for R package magick -> phytools
# python-dev-is-python3 for biopython -> superCRUNCH
# git-lfs for gittargets
# libarchive-dev for archive
# libharfbuzz-dev, libfribidi-dev for R package textshaping
# librdf0-dev for redland -> R package deposits
# cmake -> R package nanonext

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
    cd-hit \
    git-lfs \
    libarchive-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    librdf0-dev \
    libgit2-dev \
    cmake \
    wget \
    libglpk40 \
    libarchive13 \
    libzmq5 \
    liblzma-dev \
    libbz2-dev \
    libsecret-1-0 \
  && apt-get clean

########################
### python libraries ###
########################

# biopython for superCRUNCH

RUN curl https://bootstrap.pypa.io/get-pip.py | python \
  && pip install biopython

#############################
### Other custom software ###
#############################

ENV APPS_HOME=/apps
RUN mkdir $APPS_HOME
WORKDIR $APPS_HOME

### git-lfs ###
# apt-get version is too old
RUN wget curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash \
  && apt-get install git-lfs

### treePL ###
# many commits have been made since last release (v1.0)
# so checkout most recent commit (2022-04-07) 551cbde1a530adad7986c530ca1f254e3ffd42c7
ENV TPL_VERSION=551cbde1a530adad7986c530ca1f254e3ffd42c7
RUN git clone https://github.com/blackrim/treePL.git \
  && cd $APPS_HOME/treePL/ \
  && git checkout $TPL_VERSION \
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

ENV APP_NAME=gnparser
ENV GNP_VERSION=1.7.3
ENV DEST=$APPS_HOME/$APP_NAME/$GNP_VERSION
RUN wget https://github.com/gnames/gnparser/releases/download/v$GNP_VERSION/gnparser-v$GNP_VERSION-linux.tar.gz \
  && tar xf $APP_NAME-v$GNP_VERSION-linux.tar.gz \
  && rm $APP_NAME-v$GNP_VERSION-linux.tar.gz \
  && mv "$APP_NAME" /usr/local/bin/

### IQ Tree v2 ###
ENV APP_NAME=iqtree
ENV IQT_VERSION=2.2.2.7
RUN wget https://github.com/iqtree/iqtree2/releases/download/v$IQT_VERSION/iqtree-$IQT_VERSION-Linux.tar.gz \
  && tar xf $APP_NAME-$IQT_VERSION-Linux.tar.gz \
  && rm $APP_NAME-$IQT_VERSION-Linux.tar.gz \
  && mv $APP_NAME-$IQT_VERSION-Linux/bin/iqtree2 /usr/local/bin/

### trimAL ###
ENV APP_NAME=trimal
ENV TRIMAL_VERSION=492b6d9455ec2d0b19a420bcfc4eea8adc88e3ea
RUN git clone https://github.com/scapella/$APP_NAME.git && \
	cd $APP_NAME && \
  git checkout $TRIMAL_VERSION && \
  cd source && \
	make && \
	cp trimal /usr/local/bin

### taxon-tools ###
ENV APP_NAME=taxon-tools
ENV TAXONTOOLS_VERSION=5bc8c1f7b51ed51d22773f6b16d75ccec3ae6dec
RUN git clone https://github.com/camwebb/$APP_NAME.git && \
	cd $APP_NAME && \
  git checkout $TAXONTOOLS_VERSION && \
	make check && \
	make install

### pandoc ###

ENV APP_NAME=pandoc
ENV PANDOC_VERSION=3.1.6.1
RUN wget https://github.com/jgm/pandoc/releases/download/$PANDOC_VERSION/pandoc-$PANDOC_VERSION-1-amd64.deb \
  && dpkg -i pandoc-$PANDOC_VERSION-1-amd64.deb \
  && rm pandoc-$PANDOC_VERSION-1-amd64.deb

######################
### conda packages ###
######################

# install miniconda
ENV MINICONDA_VERSION py311_23.5.2-0
ENV CONDA_DIR $HOME/Miniconda3

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

### SuperCRUNCH ###
# Needs to run in a conda env
ENV APPNAME supercrunch
ENV SC_VERSION 1.3.2
ENV ENV_PREFIX /env/$APPNAME

# - Download SuperCRUNCH (python scripts and env.yaml)
RUN wget https://github.com/dportik/SuperCRUNCH/archive/refs/tags/v$SC_VERSION.tar.gz && \
  tar -xzf v$SC_VERSION.tar.gz && \
  rm v$SC_VERSION.tar.gz

# - Create conda environment
RUN conda update --name base --channel defaults conda && \
  conda install -n base conda-libmamba-solver && \
  conda config --set solver libmamba && \
  conda env create --prefix $ENV_PREFIX --file SuperCRUNCH-$SC_VERSION/$APPNAME-conda-env.yml --force && \
  conda clean --all --yes

# - Make shell script to run conda app
# this makes it so superCRUNCH scripts can be run e.g. `supercrunch Parse_Loci.py`
RUN echo '#!/bin/bash' >> /usr/local/bin/$APPNAME && \
  echo "source $CONDA_DIR/etc/profile.d/conda.sh" >> /usr/local/bin/$APPNAME && \
  echo "conda activate /env/$APPNAME" >> /usr/local/bin/$APPNAME && \
  echo "python /apps/SuperCRUNCH-$SC_VERSION/supercrunch-scripts/\"\$@\"" >> /usr/local/bin/$APPNAME && \
  chmod 755 /usr/local/bin/$APPNAME

####################################
### Install R packages with renv ###
####################################

# Create directory for renv project library
RUN mkdir /renv

# Modify Rprofile.site so renv uses /renv for project library
RUN echo 'Sys.setenv(RENV_PATHS_LIBRARY = "/renv")' >> /usr/local/lib/R/etc/Rprofile.site

# Initialize a 'dummy' project and restore the renv library.
# Since the library path is specified as above, the library will be restored to /renv
RUN mkdir /tmp/project

COPY ./renv.lock /tmp/project

WORKDIR /tmp/project

# Restore, but don't use cache
RUN Rscript -e 'install.packages("renv"); renv::consent(provided = TRUE); renv::settings$use.cache(FALSE); renv::init(bare = TRUE); renv::restore()'

############
### Cron ###
############

# cron is used to run R/setup_gb.R automatically once per day

RUN apt-get -y install cron

# Write script to launch R/setup_gb.R from /wd/
RUN echo "#!/bin/bash" >> /home/setup_gb.sh && \
  echo "cd /wd" >> /home/setup_gb.sh && \
  echo "/usr/local/bin/Rscript /wd/R/setup_gb.R" >> /home/setup_gb.sh && \
  chmod 0644 /home/setup_gb.sh

# Create the log file to be able to run tail
RUN touch /var/log/cron.log

# Setup cron job
RUN (crontab -l ; echo "0 0 * * * bash /home/setup_gb.sh >> /var/log/cron.log 2>&1") | crontab

# To run the cron job, provide the command `cron` to `docker run`:
# docker run --rm -dt -v ${PWD}:/wd -w /wd --name setup_gb joelnitta/ftol:latest cron -f
# 
# as long as the container is up, it will run the job once per day

############################
### Set up non-root user ###
############################

ARG USERNAME=ftol
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    # [Optional] Add sudo support. Omit if you don't need to install software
    # in continer.
    && apt-get update \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# Set the default user. Omit if you want to keep the default as root.
USER $USERNAME

WORKDIR /home/