#!/bin/bash
# Run workflow and pass in current docker image tag
#
# Optional first argument: a target name to pass as tar_make(names = ...),
# for a run that stops early (e.g. `run.sh non_mono_check` builds everything
# up through the taxonomy/monophyly checkpoint and no further). With no
# argument, runs the full pipeline (unchanged default behavior).
IMAGE_TAG=joelnitta/ftol:1.9.0
TARGET_NAME="${1:-}"
if [ -n "${TARGET_NAME}" ]; then
  TAR_MAKE_CALL="targets::tar_make(names = \"${TARGET_NAME}\")"
else
  TAR_MAKE_CALL="targets::tar_make()"
fi

# All tar_make() logs live in logs/, named logs/tar_make_<YYYYMMDD_HHMMSS>.log;
# logs/tar_make_latest.log always symlinks to the most recent one.
mkdir -p logs
LOG_FILE="logs/tar_make_$(date +%Y%m%d_%H%M%S).log"
ln -sf "$(basename "$LOG_FILE")" logs/tar_make_latest.log

docker run -dt --rm \
  -v ${PWD}:/wd -w /wd \
  -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) \
  -e IMAGE_TAG=$IMAGE_TAG \
  $IMAGE_TAG bash -c "Rscript -e '${TAR_MAKE_CALL}' 2>&1 | tee '$LOG_FILE'"

echo "Logging to $LOG_FILE (logs/tar_make_latest.log)"
