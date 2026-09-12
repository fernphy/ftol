#!/bin/bash
# Frequent check-in on a main-pipeline run started by main_pipeline_cron.sh.
#
# Cheap and fast when there's nothing to do (the common case): if the state
# file says idle/awaiting_approval/failed_needs_attention, or the pipeline
# still looks actively running, this exits immediately without invoking
# `claude` at all. It only calls `claude -p` once it looks like the run has
# stopped (crashed, errored, or reached the checkpoint) -- everything from
# there (diagnosing an error, applying the ppg-taxonomy-update skill,
# relaunching, deciding when to stop and notify Joel) is Claude's job per
# .claude/skills/main-pipeline-cron/SKILL.md, not duplicated here.
#
# Invoked from jnitta's crontab under flock (see below). Full write-up:
# docs/nittalab_main_pipeline_cron.md
#
#   */30 * * * * /usr/bin/flock -n /home/jnitta/ftol/.main_pipeline_monitor.lock \
#     /home/jnitta/ftol/main_pipeline_monitor.sh \
#     >> /home/jnitta/ftol/logs/main_pipeline_monitor.log 2>&1

set -uo pipefail

FTOL_DIR=/home/jnitta/ftol
STATE_FILE="${FTOL_DIR}/.main_pipeline_state"
# If _targets/meta/progress hasn't moved in this long, assume the tar_make()
# process itself has exited (finished, errored, or crashed) rather than
# still being mid-target. Should comfortably exceed this cron's own interval
# (30 min) so a slow-but-alive target isn't mistaken for a stopped run.
STALL_WINDOW_SECS=2700

cd "${FTOL_DIR}" || exit 1

state="idle"
if [ -f "${STATE_FILE}" ]; then
  # shellcheck disable=SC1090
  source "${STATE_FILE}"
fi

if [ "${state}" != "in_flight" ]; then
  # idle / awaiting_approval / failed_needs_attention: nothing for this
  # script to do. (awaiting_approval and failed_needs_attention both mean a
  # human decision is pending -- see main_pipeline_cron.sh and the skill.)
  exit 0
fi

progress="${FTOL_DIR}/_targets/meta/progress"
if [ -f "${progress}" ]; then
  age=$(( $(date +%s) - $(stat -c %Y "${progress}") ))
  if [ "${age}" -lt "${STALL_WINDOW_SECS}" ]; then
    exit 0
  fi
fi

echo "=== main_pipeline monitor: run looks stopped, invoking claude: $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="

claude -p "$(cat <<'PROMPT'
Monitor-and-fix task (see .claude/skills/main-pipeline-cron/SKILL.md,
"Monitor-and-fix mode"): the main FTOL pipeline run tracked in
.main_pipeline_state (state=in_flight) looks like it has stopped. Follow the
skill's step-by-step instructions: read the state file, inspect the targets
store per the run-tar-make skill, and either diagnose+fix+relaunch per the
ppg-taxonomy-update skill (incrementing attempt and updating the state file),
or -- if non_mono_check succeeded, or the retry cap is reached -- notify Joel
(email + push) with a summary and set state to awaiting_approval or
failed_needs_attention respectively. Do not proceed past non_mono_check under
any circumstances; that decision is Joel's alone.
PROMPT
)"
claude_status=$?

echo "=== main_pipeline monitor end (claude exit ${claude_status}): $(date -u '+%Y-%m-%d %H:%M:%S UTC') ==="
exit 0
