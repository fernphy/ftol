#!/bin/bash
# Run workflow and pass in current docker image tag
IMAGE_TAG=joelnitta/ftol:1.9.0
docker run -dt --rm \
  -v ${PWD}:/wd -w /wd \
  -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) \
  -e IMAGE_TAG=$IMAGE_TAG \
  $IMAGE_TAG Rscript -e "targets::tar_make()"
