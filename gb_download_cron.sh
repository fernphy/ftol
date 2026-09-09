#!/bin/bash
# Recurring GenBank release check for FTOL.
#
# Streams the GenBank plant division one file at a time, filters to fern +
# outgroup accessions, and rebuilds the restez database only when NCBI
# publishes a new release. Most days there is no new release, so the run
# exits non-zero within seconds (release_check errors) and builds nothing --
# that is normal, not a failure. On any non-zero exit this script prints the
# recorded target errors so the log says why (usually "No new GenBank data
# available"; the `summary` reporter set in _targets.yaml keeps the run
# itself to a single status line instead of thousands).
#
# Invoked from jnitta's crontab under flock (see below). Full write-up:
# docs/nittalab_gb_download.md
#
#   0 0 * * * /usr/bin/flock -n /home/jnitta/ftol/.gb_download.lock \
#     /home/jnitta/ftol/gb_download_cron.sh \
#     >> /home/jnitta/ftol/logs/gb_download_cron.log 2>&1

set -uo pipefail

FTOL_DIR=/home/jnitta/ftol
STORE="${FTOL_DIR}/_targets_gb_store"
# If the pipeline's progress file was touched more recently than this, assume
# a run is already in flight (possibly a manual run in another container,
# which targets' own PID-based guard does NOT see across PID namespaces) and
# skip this tick. flock in the crontab only covers cron-vs-cron overlap.
# Generous window so a long single target (DB build, tar, taxdmp) mid-run
# isn't mistaken for idle; a hard crash blocks the next tick by at most this.
ACTIVE_WINDOW_SECS=3600

cd "${FTOL_DIR}" || exit 1

echo "=== gb_download cron start: $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="

progress="${STORE}/meta/progress"
if [ -f "${progress}" ]; then
  age=$(( $(date +%s) - $(stat -c %Y "${progress}") ))
  if [ "${age}" -lt "${ACTIVE_WINDOW_SECS}" ]; then
    echo "gb_download: ${progress} modified ${age}s ago (< ${ACTIVE_WINDOW_SECS}s):" \
         "a pipeline run appears active. Skipping this tick."
    echo "=== gb_download cron end (exit 0, skipped): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
    exit 0
  fi
fi

docker_run() {
  docker run --rm \
    -v "${FTOL_DIR}":/wd -w /wd \
    -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" \
    -e TAR_PROJECT=gb_download \
    -v /mnt/jnitta/project_data/ftol_genbank_raw:/archive \
    -e GB_DL_ARCHIVE_DIR=/archive \
    joelnitta/ftol:latest \
    Rscript -e "$1"
}

docker_run 'targets::tar_make()'
status=$?

if [ "${status}" -ne 0 ]; then
  echo "--- recorded target errors (exit ${status}) ---"
  docker_run 'e <- targets::tar_meta(fields = "error", complete_only = TRUE); print(e[, c("name", "error")], n = Inf)'
fi

echo "=== gb_download cron end (exit ${status}): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
exit "${status}"
