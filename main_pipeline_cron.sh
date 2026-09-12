#!/bin/bash
# Recurring "should we start a new main FTOL analysis run" check.
#
# Unlike gb_download_cron.sh (which does everything itself), this script only
# handles the parts that are cheap and fully deterministic -- checking a new
# GenBank DB exists, PPG isn't mid-curation (main ahead of its latest release
# tag), and the previous FTOL release has actually shipped (ftol_data tagged)
# -- and launches the pipeline itself with a plain `run.sh non_mono_check`.
# Everything past that (watching the run, fixing taxonomy errors, deciding
# when to hand off to Joel) is main_pipeline_monitor.sh's job, which is why
# this script never invokes `claude`. It only does so here to send a
# notification when a condition above stops it from proceeding -- see
# notify() below.
#
# Invoked from jnitta's crontab under flock (see below). Full write-up:
# docs/nittalab_main_pipeline_cron.md
#
#   0 6 * * * /usr/bin/flock -n /home/jnitta/ftol/.main_pipeline.lock \
#     /home/jnitta/ftol/main_pipeline_cron.sh \
#     >> /home/jnitta/ftol/logs/main_pipeline_cron.log 2>&1

set -uo pipefail

FTOL_DIR=/home/jnitta/ftol
STATE_FILE="${FTOL_DIR}/.main_pipeline_state"
GB_RELEASE_FILE="${FTOL_DIR}/_targets/user/data_raw/restez/gb_release.txt"
LAST_GB_RELEASE_FILE="${FTOL_DIR}/.last_gb_release_started"
TARGETS_R="${FTOL_DIR}/_targets.R"
PPG_REPO="https://github.com/pteridogroup/ppg"
# Same reasoning as gb_download_cron.sh's ACTIVE_WINDOW_SECS, scaled up: this
# pipeline runs for days (not 1-2), and `targets`' own PID-based "already
# running" guard doesn't work across containers, so this is a backstop for a
# run started outside this automation (e.g. by hand, as documented in
# docs/nittalab_main_pipeline_cron.md -- disable this cron first). Generous
# window so a long single target mid-run isn't mistaken for idle.
ACTIVE_WINDOW_SECS=86400

cd "${FTOL_DIR}" || exit 1

echo "=== main_pipeline cron start: $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="

# notify SUBJECT BODY -- hands off to `claude -p` purely to send both
# notification channels (email via the existing send_gb_email() helper, plus
# a Claude-app push notification) with the given content, then exits. See
# .claude/skills/main-pipeline-cron/SKILL.md for the exact prompt contract.
notify() {
  local subject="$1"
  local body="$2"
  claude -p "$(cat <<PROMPT
Notify-only task (see .claude/skills/main-pipeline-cron/SKILL.md, "Notify-only
mode"): send a notification with subject "${subject}" and this body, via both
email (send_gb_email(), sourcing R/setup_gb_functions.R) and a Claude-app push
notification. Do not do anything else -- no diagnosis, no file edits, no git.

${body}
PROMPT
)" || echo "main_pipeline: notify() via claude -p failed (exit $?); see subject/body above in this log"
}

# --- Guard: is this checkout actually on main? ---
# This directory is shared with interactive use (by Joel or a Claude session)
# on this same host -- someone investigating something on a feature branch
# and forgetting to switch back is exactly the kind of accident that would
# otherwise make this script silently launch a run from the wrong code, or
# make the auto-bump-and-commit step below land a commit on the wrong branch.
current_branch="$(git -C "${FTOL_DIR}" rev-parse --abbrev-ref HEAD)"
if [ "${current_branch}" != "main" ]; then
  body="This checkout is on '${current_branch}', not main. Someone likely left"
  body="${body} it there after interactive work. Not touching it automatically --"
  body="${body} switch back to main by hand, then this will resume on its next tick."
  notify "FTOL cron: checkout is not on main" "${body}"
  echo "=== main_pipeline cron end (exit 0, skipped): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi

# --- Guard: is a run already active or awaiting a decision? ---
state="idle"
if [ -f "${STATE_FILE}" ]; then
  # shellcheck disable=SC1090
  source "${STATE_FILE}"
fi
if [ "${state}" != "idle" ]; then
  echo "main_pipeline: state=${state} (not idle). Skipping this tick."
  echo "=== main_pipeline cron end (exit 0, skipped): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi

progress="${FTOL_DIR}/_targets/meta/progress"
if [ -f "${progress}" ]; then
  age=$(( $(date +%s) - $(stat -c %Y "${progress}") ))
  if [ "${age}" -lt "${ACTIVE_WINDOW_SECS}" ]; then
    echo "main_pipeline: ${progress} modified ${age}s ago (< ${ACTIVE_WINDOW_SECS}s)" \
         "but state file says idle -- a run looks active outside this" \
         "automation's tracking (e.g. started by hand). Skipping this tick."
    echo "=== main_pipeline cron end (exit 0, skipped): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
    exit 0
  fi
fi

# --- Condition 1: is there a new local GenBank DB? ---
if [ ! -f "${GB_RELEASE_FILE}" ]; then
  echo "main_pipeline: ${GB_RELEASE_FILE} not found yet. Skipping this tick."
  echo "=== main_pipeline cron end (exit 0, skipped): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi
current_gb_release="$(cat "${GB_RELEASE_FILE}")"
last_gb_release="$(cat "${LAST_GB_RELEASE_FILE}" 2>/dev/null || echo "")"
if [ "${current_gb_release}" == "${last_gb_release}" ]; then
  echo "main_pipeline: no new GenBank release since last run (${current_gb_release})." \
       "Skipping this tick."
  echo "=== main_pipeline cron end (exit 0, skipped): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi

# --- Condition 2: is PPG's latest release trustworthy (not behind main)? ---
ppg_latest_tag="$(gh release list -R pteridogroup/ppg --limit 1 | cut -f3)"
if [ -z "${ppg_latest_tag}" ]; then
  notify "FTOL cron: couldn't check PPG release" \
    "gh release list -R pteridogroup/ppg returned nothing. Skipping today's run."
  echo "=== main_pipeline cron end (exit 0, skipped): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi
ppg_latest_sha="$(git ls-remote --tags "${PPG_REPO}" "refs/tags/${ppg_latest_tag}^{}" | cut -f1)"
if [ -z "${ppg_latest_sha}" ]; then
  # lightweight tags (no ^{} peel entry) -- fall back to the tag ref itself
  ppg_latest_sha="$(git ls-remote --tags "${PPG_REPO}" "refs/tags/${ppg_latest_tag}" | cut -f1)"
fi
ppg_main_sha="$(git ls-remote "${PPG_REPO}" main | cut -f1)"
if [ "${ppg_latest_sha}" != "${ppg_main_sha}" ]; then
  body="PPG's latest release (${ppg_latest_tag}, ${ppg_latest_sha}) does not match"
  body="${body} main (${ppg_main_sha}). main has curation ahead of the last release --"
  body="${body} please cut a new PPG release before the next automated run. Skipping"
  body="${body} today's run."
  notify "FTOL cron: PPG has unreleased changes" "${body}"
  echo "=== main_pipeline cron end (exit 0, skipped): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi

# --- Condition 3: has the previous FTOL release fully shipped? (MVP check) ---
# MVP only checks ftol_data (git tag + code= provenance); see
# docs/nittalab_main_pipeline_cron.md for the fuller checklist this is meant
# to grow into once the release process itself is automated.
( cd "${FTOL_DIR}/ftol_data" && git fetch -q origin main )
ftol_data_head="$(git -C "${FTOL_DIR}/ftol_data" rev-parse HEAD)"
ftol_data_latest_tag="$(git -C "${FTOL_DIR}/ftol_data" describe --tags --abbrev=0 2>/dev/null || echo "")"
ftol_data_tag_sha="$(git -C "${FTOL_DIR}/ftol_data" rev-list -n 1 "${ftol_data_latest_tag}" 2>/dev/null || echo "")"
if [ -z "${ftol_data_latest_tag}" ] || [ "${ftol_data_head}" != "${ftol_data_tag_sha}" ]; then
  body="ftol_data's HEAD (${ftol_data_head}) is not the same commit as its"
  body="${body} latest tag (${ftol_data_latest_tag:-none}) -- the previous release"
  body="${body} looks mid-flight (see docs/updating.md steps 9 onward). Skipping"
  body="${body} today's run."
  notify "FTOL cron: previous release not fully synced" "${body}"
  echo "=== main_pipeline cron end (exit 0, skipped): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
  exit 0
fi

# --- All clear: bump PPG version in _targets.R if needed, commit, launch ---
current_ppg_ver="$(grep -oP 'load_ppg\(ver = "\K[^"]+' "${TARGETS_R}")"
ppg_target_ver="${ppg_latest_tag#v}"
if [ "${current_ppg_ver}" != "${ppg_target_ver}" ]; then
  # Refuse to touch _targets.R if it already has uncommitted changes --
  # those are someone's in-progress interactive work, not ours to sweep into
  # an automated "Bump PPG" commit alongside our own one-line edit.
  if [ -n "$(git -C "${FTOL_DIR}" status --porcelain -- _targets.R)" ]; then
    body="Wanted to bump load_ppg(ver) from ${current_ppg_ver} to ${ppg_target_ver},"
    body="${body} but _targets.R already has uncommitted changes -- looks like"
    body="${body} in-progress interactive work. Not touching it. Skipping today's run."
    notify "FTOL cron: _targets.R has uncommitted changes" "${body}"
    echo "=== main_pipeline cron end (exit 0, skipped): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
    exit 0
  fi
  sed -i "s/load_ppg(ver = \"${current_ppg_ver}\")/load_ppg(ver = \"${ppg_target_ver}\")/" "${TARGETS_R}"
  git -C "${FTOL_DIR}" add _targets.R
  git -C "${FTOL_DIR}" commit -q -m "Bump PPG to ${ppg_latest_tag}

Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>"
  echo "main_pipeline: bumped ppg_raw ${current_ppg_ver} -> ${ppg_target_ver} and committed."
fi

echo "${current_gb_release}" > "${LAST_GB_RELEASE_FILE}"

# run.sh launches a detached (`docker run -dt`) container and returns
# immediately -- it creates its own logs/tar_make_<timestamp>.log and updates
# logs/tar_make_latest.log itself (see run.sh), so there's nothing to
# redirect or background here. Record that symlink, not a path we invent, so
# the monitor script always looks at whatever run.sh actually just started.
echo "main_pipeline: launching targeted run (through non_mono_check) via run.sh"
bash "${FTOL_DIR}/run.sh" non_mono_check

{
  echo "state=in_flight"
  echo "gb_release=${current_gb_release}"
  echo "log_file=logs/tar_make_latest.log"
  echo "attempt=0"
  echo "started_at=$(date -u '+%Y-%m-%d %H:%M:%S UTC')"
} > "${STATE_FILE}"

echo "=== main_pipeline cron end (exit 0, launched): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
exit 0
