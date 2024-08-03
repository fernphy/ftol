#!/bin/bash
# Run workflow and pass in current docker image tag
IMAGE_TAG=joelnitta/ftol:1.7.0
docker run -dt --user root --rm -v ${PWD}:/wd -w /wd --env IMAGE_TAG=$IMAGE_TAG $IMAGE_TAG Rscript -e "targets::tar_make()"
