---
name: main-pipeline-cron
description: Handle a `claude -p` invocation from main_pipeline_cron.sh or main_pipeline_monitor.sh (the cron automation that starts and watches unattended main-FTOL-pipeline runs on the nittalab server). Covers two modes -- sending a plain notification, and the monitor-and-fix loop that diagnoses a stopped run, applies ppg-taxonomy-update/run-tar-make fixes, relaunches up to a retry cap, and always stops for Joel's explicit go-ahead before sanger_ml_tree_rep_* (the 7-10 day ML tree stage). Use when invoked non-interactively by either of those two scripts, or when picking up where an "awaiting_approval"/"failed_needs_attention" cron run left off.
user-invocable: true
---

# main-pipeline-cron — unattended pipeline runs, with a human checkpoint before the ML tree

This skill is the "what Claude does" half of a two-script cron automation
(`main_pipeline_cron.sh` daily, `main_pipeline_monitor.sh` every ~30 min — see
`docs/nittalab_main_pipeline_cron.md` for the full setup). Both scripts do their own
deterministic checks in plain bash and only invoke `claude -p` for the parts that need
judgment: sending a notification, or diagnosing/fixing/deciding what a stopped run
means. This skill is purely about *that* part — for the actual taxonomy mechanics use
`ppg-taxonomy-update`, and for the actual launch/log mechanics use `run-tar-make`; don't
duplicate either here.

**Hard rule, no exceptions: never advance the pipeline past `non_mono_check`** (i.e.
never run plain `run.sh` / untargeted `tar_make()`) as part of this skill. The next
stage, `sanger_ml_tree_rep_*`, takes 7–10 days — starting it is always Joel's explicit
decision, made after reading a summary, never something this automation decides on its
own regardless of how clean the run looks.

## State file: `.main_pipeline_state`

Plain `KEY=value` lines at the repo root (gitignored), written by
`main_pipeline_cron.sh` and updated by this skill. `source`-able from bash.

| Key | Meaning |
|---|---|
| `state` | `idle` \| `in_flight` \| `awaiting_approval` \| `failed_needs_attention` |
| `gb_release` | GenBank release number this run is for |
| `log_file` | Path to the current `tar_make()` log (normally `logs/tar_make_latest.log`) |
| `attempt` | Number of diagnose-fix-relaunch cycles used so far |
| `started_at` | UTC timestamp the run was launched |

`state` values and who's allowed to touch them:

- `idle` — set only by a human (or a session acting on Joel's explicit instruction)
  after `awaiting_approval`/`failed_needs_attention` has been reviewed and resolved.
  `main_pipeline_cron.sh` will only start a new run when it sees this.
- `in_flight` — set by `main_pipeline_cron.sh` at launch. This skill may relaunch and
  stay in this state (incrementing `attempt`), or move it to `awaiting_approval` /
  `failed_needs_attention`. Never move it back to `idle` yourself.
- `awaiting_approval` — `non_mono_check` succeeded. Nothing to fix; waiting on Joel.
- `failed_needs_attention` — retry cap reached without reaching `non_mono_check`.
  Waiting on Joel to look at it (possibly by hand, possibly by asking a Claude session
  to keep trying with a fresh approach).

**Retry cap**: default 5 attempts. If `docs/nittalab_main_pipeline_cron.md` documents a
different value for this install, use that instead.

## Notify-only mode

Invoked by `main_pipeline_cron.sh` when one of its own checks fails (PPG drift, GenBank
release lookup failure, previous release not synced). The prompt gives you an exact
subject and body. Do exactly this, nothing more:

1. Send the email:
   ```r
   source("R/setup_gb_functions.R")
   send_gb_email(subject = "<given subject>", body_html = "<given body>")
   ```
   (`send_gb_email()` is a generic subject/body sender despite its name — shared with
   the GenBank-download cron, not GenBank-specific in what it sends.)
2. Send a Claude-app push notification with the same content via the `PushNotification`
   tool.
3. Exit. Don't inspect the pipeline, don't touch git, don't touch the state file — this
   mode is for conditions `main_pipeline_cron.sh` already fully diagnosed itself.

## Monitor-and-fix mode

Invoked by `main_pipeline_monitor.sh` once `_targets/meta/progress` has gone quiet long
enough that the `tar_make()` process it launched looks like it has stopped (finished,
errored, or crashed) — but the script doesn't know *why*, only that it's worth looking.

1. **Read `.main_pipeline_state`.** Confirm `state=in_flight`; if not, something raced
   with another invocation — log a note and exit without acting.
2. **Assess what actually happened**, using the techniques in the `run-tar-make` skill
   (`tar_meta(fields = "error", complete_only = TRUE)`, `tar_progress()`,
   `tar_workspace(<target>)`) against this repo's main `_targets` store:
   - **`non_mono_check` completed with no error** → success. Go to *Checkpoint reached*
     below.
   - **Some target errored** → go to *Diagnose and fix* below.
   - **Nothing errored and `non_mono_check` hasn't completed, but the process really is
     gone** (no live `targets::tar_make` process reachable, and `_targets/meta/process`
     names a dead/zombie pid per the `run-tar-make` skill's zombie note) → treat as an
     unexplained crash: log what you found, then follow *Diagnose and fix* below (relaunch
     if under the cap) since a resumed `tar_make()` will just pick up where it left off.
3. **Diagnose and fix** (only reached for an actual error or unexplained crash):
   - If `attempt >= cap`: skip straight to *Cap reached* below — do not attempt another
     fix.
   - Otherwise, diagnose the error using the `ppg-taxonomy-update` skill if it's a
     taxonomy/name-resolution failure (the large majority of cases this automation
     exists for). If it's clearly *not* a taxonomy problem (e.g. an infra issue like the
     historical `mpcheck_monophy` hang, already fixed, or something novel), use your
     best judgment and existing FTOL context to fix it the same way you would
     interactively — but if you're not confident a fix is correct (e.g. it would require
     a taxonomic judgment call only Joel can make, mirroring the "keep original name vs.
     comb. ined." decision in `ppg-taxonomy-update`), stop and treat it like a cap-hit:
     notify with what you found and set `failed_needs_attention`, rather than guessing.
   - Apply the fix, verify it standalone exactly as `ppg-taxonomy-update` describes
     (never skip straight to relaunching on a hunch), commit it (same attribution
     convention as every other commit this project uses:
     `Co-Authored-By: Claude Sonnet 5 <noreply@anthropic.com>`), then relaunch:
     `bash run.sh non_mono_check`. Update the state file: increment `attempt`, keep
     `log_file=logs/tar_make_latest.log` (run.sh always points that symlink at the new
     run), leave `state=in_flight`. Exit — the *next* `main_pipeline_monitor.sh` tick
     will pick up from here once this new run also goes quiet.
4. **Checkpoint reached** (`non_mono_check` succeeded): set `state=awaiting_approval`.
   Notify (email via `send_gb_email()` + push via `PushNotification`) with a summary —
   see *What the summary should contain* below. Exit. Do not relaunch, do not run the
   full pipeline.
5. **Cap reached** (attempted `attempt >= cap` fixes without reaching `non_mono_check`):
   set `state=failed_needs_attention`. Notify (both channels) with a summary of every
   attempt made and the current/last error. Exit.

### What the summary should contain

Whether it's a success (*Checkpoint reached*) or a cap-hit (*Cap reached*), Joel is
reading this without having watched any of it happen, so be concrete, not just "done" or
"failed":

- GenBank release this run is for, and how many attempts it took (or that it succeeded
  cleanly on the first try).
- For each fix applied: which target errored, what the root cause was, what changed
  (files, and a one-line rationale — not a full diff), and whether an FTOL/PPG issue was
  filed for it (per `ppg-taxonomy-update`'s issue-filing guidance).
- Current state (`awaiting_approval` or `failed_needs_attention`) and what happens next:
  for `awaiting_approval`, that running `run.sh` (no target argument, the full pipeline)
  is all that's needed once Joel says go; for `failed_needs_attention`, what's still
  broken and what you'd try next if asked to keep going.

## Resuming after Joel approves

Not something this skill does automatically — Joel reviews the summary and tells
*whichever* Claude session he's talking to (interactive, this one, or another cron
invocation he triggers by hand) to proceed. That session should: confirm
`state=awaiting_approval`, set `state=idle` in `.main_pipeline_state`, and run `run.sh`
(no argument — the full pipeline) per the `run-tar-make` skill. From there
`main_pipeline_monitor.sh` has nothing further to do (`sanger_ml_tree_rep_*` and
everything after it isn't taxonomy-error-prone the way name resolution is, and isn't
this skill's concern) — Joel checks in on that run the normal way (`run-tar-make`'s
monitoring section), on his own timeline, since it's the 7–10 day step and there's no
further automated checkpoint past it in this design.
