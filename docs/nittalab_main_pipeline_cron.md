# Setting up automated main-pipeline runs on the nittalab server

This is a guide for Claude (or a human) implementing recurring, automated main-FTOL-
pipeline runs on the nittalab server, following the same host-crontab pattern as
GenBank downloads (`docs/nittalab_gb_download.md`). Read that doc first if you haven't
— this one assumes its conventions (flock, staleness guards, `docker run` per
invocation, enable/disable via `crontab -l | sed ... | crontab -`).

**Do not enable this until you've confirmed no `tar_make()` is currently running
against this repo's `_targets` store** (check `_targets/meta/progress`'s mtime, and ask
Joel if in doubt — as of when this was written, a manual run was in progress). Installing
the crontab entries below is safe to do at any time (they're inert until their first
scheduled tick); just don't let the first tick fire while something else is using the
store.

**None of this exists on `main` yet.** It was developed on `feat/main-pipeline-cron`
(PR: check for one titled "Add cron automation for the main FTOL pipeline" if it's not
linked from wherever you found this doc) specifically so the in-progress run above
wouldn't be touched by anything in this checkout while it was written. Nothing here —
scripts, this doc, the `main-pipeline-cron` skill — can be installed or relied on until
that PR is merged to `main`, since the host crontab entries below point at scripts that
only exist post-merge. Also: this checkout is shared between interactive use (Joel or a
Claude session working in a devcontainer bind-mounted to the same path) and whatever
cron runs here — if you ever find this checkout on a branch other than `main` outside of
active feature work, that's a mistake to fix (`git checkout main`), not a new normal;
`main_pipeline_cron.sh` itself refuses to proceed (and notifies) if it isn't.

## Why this is two scripts, not one

GenBank downloads are a single deterministic script because there's nothing to decide —
either NCBI has a new release or it doesn't, and the pipeline either builds cleanly or
it doesn't. The main pipeline is different: it reliably hits taxonomy/name-resolution
errors that need a taxonomic judgment call to fix (see `ppg-taxonomy-update`), and its
last stage (`sanger_ml_tree_rep_*`) takes 7–10 days, so it needs a deliberate human
checkpoint rather than running end-to-end unattended. That splits into:

- **`main_pipeline_cron.sh`** (daily) — purely deterministic: is there a new GenBank
  release, is PPG in a releasable state, has the previous FTOL release shipped. All
  three are plain `git`/`gh`/file checks, no judgment needed, so this script never
  invokes `claude`. If everything's clear it launches a *targeted* run
  (`run.sh non_mono_check` — stops at the taxonomy/monophyly checkpoint, well short of
  the ML tree stage) and exits.
- **`main_pipeline_monitor.sh`** (every ~30 min) — cheap to run when idle (exits
  immediately without touching `claude` if nothing's in flight or nothing's changed).
  Once a launched run looks like it's stopped, it hands off to `claude -p`, which
  diagnoses what happened, fixes it via `ppg-taxonomy-update`/`run-tar-make`, and
  relaunches — up to a retry cap — or stops and notifies Joel. See
  `.claude/skills/main-pipeline-cron/SKILL.md` for exactly what Claude does in that
  invocation; this doc is just the cron/host plumbing around it.

Both scripts are unconditionally committed at the repo root, same as `gb_download_cron.sh`
and `run.sh`.

## One-time setup

1. Everything `docs/nittalab_gb_download.md`'s one-time setup already covers (stable
   repo checkout path, `.secrets/` gmailr creds, the `joelnitta/ftol` image built) is
   shared — nothing new needed there.
2. **`claude` CLI installed and authenticated** on the host (not just inside a
   devcontainer) as the account that should own these automated runs, and confirmed
   runnable non-interactively: `claude -p "say hi"` should print a response and exit
   with no prompts. Cron invocations can't answer an interactive permission prompt, so
   confirm whatever permission mode is needed for it to edit files/run bash/use git
   without pausing (check `claude --help` for the current flag — this changes between
   CLI versions, so don't assume the exact name; `main_pipeline_monitor.sh`'s `claude
   -p ...` calls in this repo don't currently pass one, add it there once confirmed).
3. **`gh` CLI authenticated** on the host (used for the PPG release check and, if you
   ever extend the sync check past the ftol_data MVP, GitHub API lookups on other
   repos). Check with `gh auth status` — don't assume a working devcontainer setup
   carries over to the bare host or a fresh container image. If it's not authenticated,
   `gh auth login --hostname github.com --web` gives a device code you approve at
   github.com/login/device from any browser (answer "No" if it offers to generate a new
   SSH key — that's separate from `gh`'s own auth and not needed if git push already
   works). Also confirm `git push`/`git ls-remote` actually work — this repo's remote is
   SSH (`git@github.com:...`), which needs `openssh-client` installed and GitHub's host
   key in `~/.ssh/known_hosts` (`ssh-keyscan -H github.com >> ~/.ssh/known_hosts`);
   both were missing the first time this was tried in the devcontainer this was
   developed in, so don't assume either is already in place on a fresh host/container.
4. **A local `ftol_data` clone**, fetchable, at `<repo>/ftol_data` (already the case in
   the devcontainer this was developed in; confirm it also exists on the host at the
   same relative path — `main_pipeline_cron.sh` assumes `${FTOL_DIR}/ftol_data`).
5. **Push notifications**: confirm the `PushNotification` tool actually reaches Joel's
   Claude app from a `claude -p` invocation launched by cron (not just from an
   interactive session) — this is untested as of writing. If it doesn't work in this
   context, `.claude/skills/main-pipeline-cron/SKILL.md`'s "Notify-only mode" and
   "Monitor-and-fix mode" both fall back gracefully to email-only; no code change is
   needed, just note it in this doc once confirmed either way.

## Host crontab entries

```cron
0 6 * * * /usr/bin/flock -n /home/jnitta/ftol/.main_pipeline.lock /home/jnitta/ftol/main_pipeline_cron.sh >> /home/jnitta/ftol/logs/main_pipeline_cron.log 2>&1
*/30 * * * * /usr/bin/flock -n /home/jnitta/ftol/.main_pipeline_monitor.lock /home/jnitta/ftol/main_pipeline_monitor.sh >> /home/jnitta/ftol/logs/main_pipeline_monitor.log 2>&1
```

- The daily check runs at 06:00 — a few hours after gb_download's midnight slot, so a
  same-day new GenBank release has time to finish downloading (`gb_release.txt` is
  deliberately the very last file that pipeline touches) before this one looks for it.
  On a day with nothing new, this just means a one-day delay, not a missed run.
- The monitor runs every 30 minutes. It's cheap when idle (one file read, maybe one
  `stat`) — no `claude` invocation unless there's actually something to look at.
- Each has its own lock file (`.main_pipeline.lock`, `.main_pipeline_monitor.lock`) so
  they never block each other; the *pipeline itself* not double-starting is guarded by
  `.main_pipeline_state` (see the skill) plus the same `_targets/meta/progress`
  staleness backstop `gb_download_cron.sh` uses for cross-container overlap.
- `logs/` must exist and be writable (gitignored, same as the gb_download logs).

## Enabling and disabling both jobs

Same non-interactive `crontab -l | sed ... | crontab -` pattern as
`docs/nittalab_gb_download.md`. Disable **both** lines before running the pipeline by
hand, and whenever the server is being worked on:

```bash
crontab -l | sed \
  -e 's|^0 6 \* \* \*|#DISABLED# 0 6 * * *|' \
  -e 's|^\*/30 \* \* \* \*|#DISABLED# */30 * * * *|' \
  | crontab -
```

Re-enable:

```bash
crontab -l | sed 's|^#DISABLED# ||' | crontab -
```

Confirm with `crontab -l`.

## What to expect

- Most days: `main_pipeline_cron.sh` finds no new GenBank release and exits
  immediately, logging why. `main_pipeline_monitor.sh` finds `state=idle` and exits
  immediately, logging nothing (it's silent on the no-op path — check
  `.main_pipeline_state` directly if you want to confirm it's running at all).
- When a new release triggers a run: expect it to reach `non_mono_check` within hours
  (per this session's experience: name resolution is minutes, the compute phase
  including monophyly checks and Sanger/plastome tree-building is a few hours), *not*
  days — the multi-day step is `sanger_ml_tree_rep_*`, which this automation never
  starts on its own. You'll get an email + push notification either when it reaches
  that checkpoint cleanly, or when the retry cap (default 5 attempts — see the skill)
  is hit without getting there. Either way, check `.main_pipeline_state` and the
  summary in the notification, then tell a Claude session to proceed (or to keep
  investigating) when you're ready — see the skill's "Resuming after Joel approves".
- If PPG has unreleased changes (main ahead of its latest tag) or the previous FTOL
  release hasn't fully shipped (MVP check: `ftol_data` not tagged at its current HEAD),
  you'll get a notification and the day is skipped — no run started, `state` stays
  `idle`, tomorrow's tick tries again.
- If a run gets interrupted (host reboot, container killed) with `state=in_flight` and
  no error recorded, the monitor treats that as an unexplained crash on its next tick
  and relaunches (targets' own caching means completed targets aren't redone) — same
  resumability property as gb_download.

## Debugging

- `.main_pipeline_state` (repo root, gitignored) is the single source of truth for what
  this automation thinks is happening — read it directly rather than inferring from
  logs.
- `logs/main_pipeline_cron.log` / `logs/main_pipeline_monitor.log` are these two
  scripts' own output (bracketed start/end markers, same style as
  `logs/gb_download_cron.log`). The actual `tar_make()` output is wherever
  `.main_pipeline_state`'s `log_file` points (normally `logs/tar_make_latest.log` — see
  `run-tar-make`).
- To inspect the pipeline directly, use the techniques in `run-tar-make` and
  `ppg-taxonomy-update` against this repo's main `_targets` store — nothing about the
  cron wrapper changes how the store itself is inspected.
- To force a retry from `failed_needs_attention` without waiting for a code fix (e.g.
  you fixed something by hand), reset `.main_pipeline_state`'s `attempt` to `0` and
  `state` to `in_flight`; the next monitor tick will pick it up as if it had just gone
  quiet.
