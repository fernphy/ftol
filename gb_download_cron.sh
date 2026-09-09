#!/bin/bash
# Recurring GenBank release check for FTOL.
#
# Streams the GenBank plant division one file at a time, filters to fern +
# outgroup accessions, and rebuilds the restez database only when NCBI
# publishes a new release. Most days this exits non-zero almost immediately
# with "No new GenBank data available" -- that is normal, not a failure.
#
# Invoked from jnitta's crontab under flock (see below). Full write-up:
# docs/nittalab_gb_download.md
#
#   0 0 * * * /usr/bin/flock -n /home/jnitta/ftol/.gb_download.lock \
#     /home/jnitta/ftol/gb_download_cron.sh \
#     >> /home/jnitta/ftol/logs/gb_download_cron.log 2>&1

set -uo pipefail

cd /home/jnitta/ftol || exit 1

echo "=== gb_download cron start: $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="

docker run --rm \
  -v /home/jnitta/ftol:/wd -w /wd \
  -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" \
  -e TAR_PROJECT=gb_download \
  -v /mnt/jnitta/project_data/ftol_genbank_raw:/archive \
  -e GB_DL_ARCHIVE_DIR=/archive \
  joelnitta/ftol:latest \
  Rscript -e 'targets::tar_make()'
status=$?

echo "=== gb_download cron end (exit ${status}): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
exit "${status}"
