#!/bin/bash
# Recurring driver for the FTOL main pipeline, chained after
# gb_download_cron.sh has (maybe) pulled a new GenBank release.
#
# Detects new data by comparing the current_release target (the release
# gb_download has actually pulled) recorded in _targets_gb_store against the
# gb_release target in the main _targets store; if they
# diverge, runs run.sh (never a hand-rolled docker run -- see the
# run-tar-make skill), waits for it to finish, then runs
# R/publish_figshare.R and R/snapshot_ftol_data.R headlessly and emails
# joelnitta@gmail.com either way. Stops there -- version bumps/releases in
# ftol_data/ftolr/ftol_vis/ftol_shiny/the website always need a human (see
# the release-ftol skill).
#
# Invoked from jnitta's crontab under flock (see below). Full write-up:
# docs/nittalab_main_pipeline_cron.md
#
#   0 3 * * * /usr/bin/flock -n /home/jnitta/ftol/.main_pipeline.lock \
#     /home/jnitta/ftol/main_pipeline_cron.sh \
#     >> /home/jnitta/ftol/logs/main_pipeline_cron.log 2>&1

set -uo pipefail

FTOL_DIR=/home/jnitta/ftol
GITCONFIG_HOST_FILE=/home/jnitta/.gitconfig
# Run containers as the host user (entrypoint.sh remaps via HOST_UID/GID) so
# nothing they write into the bind-mounted repos ends up root-owned.
HOST_USER_ARGS="-e HOST_UID=$(id -u) -e HOST_GID=$(id -g)"
ACTIVE_WINDOW_SECS=3600

cd "${FTOL_DIR}" || exit 1

echo "=== main_pipeline cron start: $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="

# Staleness guard, same rationale as gb_download_cron.sh: targets' PID-based
# "already running" guard doesn't cross container/PID-namespace boundaries.
progress="${FTOL_DIR}/_targets/meta/progress"
if [ -f "${progress}" ]; then
  age=$(( $(date +%s) - $(stat -c %Y "${progress}") ))
  if [ "${age}" -lt "${ACTIVE_WINDOW_SECS}" ]; then
    echo "main_pipeline: ${progress} modified ${age}s ago (< ${ACTIVE_WINDOW_SECS}s):" \
         "a pipeline run appears active. Skipping this tick."
    echo "=== main_pipeline cron end (exit 0, skipped): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
    exit 0
  fi
fi

docker_read() {
  docker run --rm ${HOST_USER_ARGS} -v "${FTOL_DIR}":/wd -w /wd \
    joelnitta/ftol:latest Rscript -e "$1"
}

# The download store's current_release and the main store's gb_release both
# hold the GenBank release number; they diverge exactly
# when gb_download has pulled a release the main pipeline hasn't processed
# yet. No separate marker file needed.
gb_release_downloaded=$(docker_read \
  'cat(targets::tar_read(current_release, store = "_targets_gb_store"))' 2>/dev/null)
gb_release_processed=$(docker_read \
  'cat(targets::tar_read(gb_release, store = "_targets"))' 2>/dev/null)

if [ -z "${gb_release_downloaded}" ]; then
  echo "main_pipeline: could not read current_release from _targets_gb_store (has" \
       "gb_download_cron.sh completed a run yet?). Skipping."
  echo "=== main_pipeline cron end (exit 0, no gb_download data): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi

if [ "${gb_release_downloaded}" = "${gb_release_processed}" ]; then
  echo "main_pipeline: gb_release unchanged (${gb_release_processed}). Nothing to do."
  echo "=== main_pipeline cron end (exit 0, no new data): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi

echo "main_pipeline: new GenBank release ${gb_release_downloaded} (main pipeline last" \
     "processed ${gb_release_processed:-none}). Running run.sh."

# IMAGE_TAG is hard-coded in run.sh and bumped by hand per release; read it
# from there rather than duplicating/hardcoding it here.
image_tag=$(grep -oP '(?<=^IMAGE_TAG=)\S+' run.sh | head -1)
bash run.sh

# run.sh launches a detached, --rm container; wait for it to exit.
sleep 10
while docker ps --filter "ancestor=${image_tag}" --format '{{.ID}}' | grep -q .; do
  sleep 60
done

error_count=$(docker_read \
  'e <- targets::tar_meta(fields = "error", complete_only = TRUE); cat(nrow(e))')

if [ "${error_count}" != "0" ]; then
  echo "main_pipeline: tar_make() reported ${error_count} error(s); see" \
       "logs/tar_make_latest.log."
  docker_read "source('R/setup_gb_functions.R'); send_release_ready_email(status = 'pipeline_failed')"
  echo "=== main_pipeline cron end (exit 1, pipeline failed): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 1
fi

echo "main_pipeline: tar_make() finished without error. Running snapshot + FigShare publish."

# The gitconfig is only for the commit author; the container has no gh and no
# credentials, so it can't push -- the host does that below.
docker run --rm ${HOST_USER_ARGS} \
  -v "${FTOL_DIR}":/wd -w /wd \
  -v "${GITCONFIG_HOST_FILE}":/etc/gitconfig_persisted:ro \
  -e GIT_CONFIG_GLOBAL=/etc/gitconfig_persisted \
  joelnitta/ftol:latest \
  Rscript -e "source('R/publish_figshare.R'); source('R/snapshot_ftol_data.R')"
publish_status=$?

if [ "${publish_status}" -eq 0 ]; then
  # Push the snapshot commit from the host, which has the gh credentials
  # (no-op "Everything up-to-date" if the snapshot had nothing to commit).
  git -C "${FTOL_DIR}/ftol_data" push origin main
  publish_status=$?
fi

if [ "${publish_status}" -eq 0 ]; then
  docker_read "source('R/setup_gb_functions.R'); send_release_ready_email(status = 'ready', gb_release = ${gb_release_downloaded})"
  echo "=== main_pipeline cron end (exit 0, release ready): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
else
  echo "main_pipeline: publish_figshare.R / snapshot_ftol_data.R failed (exit ${publish_status})."
  docker_read "source('R/setup_gb_functions.R'); send_release_ready_email(status = 'publish_failed', gb_release = ${gb_release_downloaded})"
  echo "=== main_pipeline cron end (exit 1, publish failed): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 1
fi
