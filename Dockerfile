FROM rocker/rstudio:4.0.0

ARG DEBIAN_FRONTEND=noninteractive

############################
### Install APT packages ###
############################

# gcc through libtool for treePL
# cmake, libeigen3-dev for IQTREE
# zlib1g-dev for R package XVector
# libxml2-dev for R package XML
# libudunits2-dev for R package units
# libgdal-dev for R package sf
# libzmq3-dev for R package rzmq -> clustermq
# libmagick++-dev for R package magick -> phytools
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
  && apt-get clean

####################################
### Install R packages with renv ###
####################################

# Create directory for renv project library
RUN mkdir renv

# Modify Rprofile.site so renv uses /renv for project library
RUN echo 'Sys.setenv(RENV_PATHS_LIBRARY = "/renv")' >> /usr/local/lib/R/etc/Rprofile.site

# Initialize a 'dummy' project and restore the renv library.
# Since the library path is specified as above, the library will be restored to /renv
RUN mkdir tmp/project

COPY ./renv.lock tmp/project

WORKDIR tmp/project

# Don't use cache (the symlinks won't work from Rstudio server)
RUN Rscript -e 'install.packages("renv"); renv::consent(provided = TRUE); renv::settings$use.cache(FALSE); renv::init(bare = TRUE); renv::restore()'

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
ENV VERSION=0.14.1
ENV DEST=$APPS_HOME/$APP_NAME/$VERSION
RUN wget https://gitlab.com/gogna/gnparser/uploads/7d6ed7e3b1eee0fd6c9ae51f5bf711c0/$APP_NAME-v$VERSION-linux.tar.gz \
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

WORKDIR /home/rstudio/
