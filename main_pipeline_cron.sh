#!/bin/bash
# Recurring driver for the FTOL main pipeline, chained after
# gb_download_cron.sh has (maybe) pulled a new GenBank release.
#
# Detects new data by comparing the current_release target (the release
# gb_download has actually pulled) recorded in _targets_gb_store against the
# last release this wrapper fully processed (_targets/user/
# .last_processed_gb_release, written only after pipeline + publish + push all
# succeed, so a failed run is retried); if they
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
# GenBank release that has been *fully* processed: pipeline + FigShare publish
# + snapshot + push all succeeded. Written only at the very end of a good run,
# so a run that fails partway is retried on the next tick. (The main store's
# own gb_release target can't serve as this marker: it completes in the first
# seconds of any run, so after a mid-run failure it already says "processed".)
# Lives in the gitignored _targets/user/ so it never dirties the repo. After
# processing a release by hand, update it: echo <release> > this file.
STATE_FILE="${MAIN_PIPELINE_STATE_FILE:-${FTOL_DIR}/_targets/user/.last_processed_gb_release}"
# MAIN_PIPELINE_DRY_RUN=1: decide and report, but never start a run
DRY_RUN="${MAIN_PIPELINE_DRY_RUN:-0}"

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

# Second guard: the progress file above is only touched when a target starts or
# finishes, so during a long single target (the ~10-day IQ-TREE replicates) it
# goes quiet for hours and the guard passes while a run is still going. A
# container from the run.sh image still running means a run is active, whether
# cron or a human launched it. (Two runs share one _targets store and one set
# of IQ-TREE checkpoints and corrupt each other.)
# IMAGE_TAG is hard-coded in run.sh and bumped by hand per release; read it
# from there rather than duplicating/hardcoding it here.
image_tag=$(grep -oP '(?<=^IMAGE_TAG=)\S+' run.sh | head -1)
if docker ps --filter "ancestor=${image_tag}" --format '{{.ID}}' | grep -q .; then
  echo "main_pipeline: a container from ${image_tag} is running (a pipeline" \
       "run is active). Skipping this tick."
  echo "=== main_pipeline cron end (exit 0, skipped): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi

docker_read() {
  docker run --rm ${HOST_USER_ARGS} -v "${FTOL_DIR}":/wd -w /wd \
    joelnitta/ftol:latest Rscript -e "$1"
}

# gb_download's current_release (the release it has actually pulled) vs the
# last release this wrapper fully processed (STATE_FILE).
gb_release_downloaded=$(docker_read \
  'cat(targets::tar_read(current_release, store = "_targets_gb_store"))' 2>/dev/null)
if [ -s "${STATE_FILE}" ]; then
  gb_release_processed=$(tr -d '[:space:]' < "${STATE_FILE}")
else
  # Deliberately NOT guessed from the main store: a wrong guess would start a
  # run that publishes to FigShare. Seed it once by hand with the last release
  # that was fully processed and released.
  echo "main_pipeline: ${STATE_FILE} is missing or empty. Create it with the" \
       "last fully processed GenBank release (e.g. echo 273 > that file)." \
       "Refusing to guess. Skipping."
  echo "=== main_pipeline cron end (exit 1, no state file): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 1
fi

if [ -z "${gb_release_downloaded}" ]; then
  echo "main_pipeline: could not read current_release from _targets_gb_store (has" \
       "gb_download_cron.sh completed a run yet?). Skipping."
  echo "=== main_pipeline cron end (exit 0, no gb_download data): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi

if [ "${gb_release_downloaded}" = "${gb_release_processed}" ]; then
  echo "main_pipeline: release ${gb_release_processed} already fully processed. Nothing to do."
  echo "=== main_pipeline cron end (exit 0, no new data): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi

echo "main_pipeline: new GenBank release ${gb_release_downloaded} (last fully" \
     "processed: ${gb_release_processed:-none})."
if [ "${DRY_RUN}" = "1" ]; then
  echo "main_pipeline: DRY RUN -- would run run.sh now. Stopping."
  echo "=== main_pipeline cron end (exit 0, dry run): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi
echo "main_pipeline: Running run.sh."

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
  # Only now is this release "done"; a failure anywhere above leaves the marker
  # untouched so the next tick retries (targets resumes from cache).
  printf '%s\n' "${gb_release_downloaded}" > "${STATE_FILE}"
  date_cutoff=$(docker_read 'cat(as.character(targets::tar_read(date_cutoff, store = "_targets")))' 2>/dev/null)
  docker_read "source('R/setup_gb_functions.R'); send_release_ready_email(status = 'ready', gb_release = ${gb_release_downloaded}, date_cutoff = '${date_cutoff:-unknown}')"
  echo "=== main_pipeline cron end (exit 0, release ready): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
else
  echo "main_pipeline: publish_figshare.R / snapshot_ftol_data.R failed (exit ${publish_status})."
  docker_read "source('R/setup_gb_functions.R'); send_release_ready_email(status = 'publish_failed', gb_release = ${gb_release_downloaded})"
  echo "=== main_pipeline cron end (exit 1, publish failed): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 1
fi
