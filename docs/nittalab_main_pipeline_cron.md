# Setting up automated main-pipeline runs on the nittalab server

This is a guide for Claude (or a human) enabling `main_pipeline_cron.sh` — the
cron job that chains the main FTOL pipeline (`run.sh`) and the data
snapshot/FigShare publish steps after `gb_download_cron.sh` has pulled a new
GenBank release. Companion to `docs/nittalab_gb_download.md`, which covers
the download-only job this one is chained after.

## Why this exists

`gb_download_cron.sh` (see `docs/nittalab_gb_download.md`) only runs the
`gb_download` project (`_targets_gb.R`). Nothing currently triggers the
`main` project (`_targets.R`, i.e. `run.sh`'s target) automatically — a human
has to notice a new release and run `run.sh` by hand, then separately run
`R/publish_figshare.R` and `R/snapshot_ftol_data.R`. `main_pipeline_cron.sh`
closes that gap, stopping right before the point a human decision is
actually needed (bumping/releasing versions in `ftol_data`, `ftolr`,
`ftol_vis`, `ftol_shiny`, or the website — see `docs/updating.md` /
the `release-ftol` skill).

## One-time setup

1. Confirm the **host** can push to GitHub non-interactively, since the host
   (not the container) pushes `ftol_data` after the snapshot: `gh auth
   status` should show a login, and `~/.gitconfig` should have the
   `credential.helper = !/usr/bin/gh auth git-credential` entries (from
   `gh auth setup-git`, with git protocol `https`) so `git push` needs no SSH
   agent or live terminal. Verify with `git -C ftol_data push --dry-run
   origin main` from a plain shell (the cron environment uses the same
   `~/.gitconfig`; a VS Code terminal can mask gaps via its own SSH agent
   forwarding).

   The release container gets none of these credentials: it has no `gh`, and
   mounts `~/.gitconfig` read-only at `/etc/gitconfig_persisted`
   (`GIT_CONFIG_GLOBAL`) only so commits carry the right author. That file
   should also contain `[safe] directory = *` (a fresh container has no
   other safe.directory entries for the mounted repos). The persisted
   `/home/jnitta/.gh_config` (`GH_CONFIG_DIR`) is only used by the dev
   container (see `docker-compose.yml`), not by this cron job.
2. Confirm `FIGSHARE_TOKEN` is present in `.Renviron` (already required for
   `R/publish_figshare.R`, already used by the pre-existing
   `upload_to_figshare()` helper).
3. Confirm `.secrets/` (gmailr OAuth credentials) is present — same
   credentials `gb_download_cron.sh` already uses for its start/finish
   emails, reused here by `send_release_ready_email()`
   (`R/setup_gb_functions.R`).
4. **Before enabling the crontab entry below, check for a job nobody
   remembered**: this container can't read `root`'s or another user's
   crontab (`sudo crontab -l` / `sudo ls /var/spool/cron/crontabs` needs a
   human on the host). Confirmed so far: `jnitta`'s crontab only runs
   `gb_download_cron.sh`; the `main` project has no visible cron/timer.

## Host crontab entry

Scheduled a few hours after the existing `gb_download_cron.sh` tick (daily
`0 0 * * *` JST) to let it finish first — most days it exits in ~2.5 min, but
budget for the rare multi-day GenBank-release-processing run:

```cron
0 3 * * * /usr/bin/flock -n /home/jnitta/ftol/.main_pipeline.lock /home/jnitta/ftol/main_pipeline_cron.sh >> /home/jnitta/ftol/logs/main_pipeline_cron.log 2>&1
```

`logs/` must exist and be writable by `jnitta` (gitignored, shared with
`gb_download_cron.sh`'s log). `flock` creates `.main_pipeline.lock` on first
run (gitignored) — a separate lock file from `gb_download_cron.sh`'s, since
these are two independent jobs.

[`main_pipeline_cron.sh`](../main_pipeline_cron.sh) (repo root, committed
alongside `run.sh`/`gb_download_cron.sh`):

1. Skips the tick (exit 0) if `_targets/meta/progress` was touched in the
   last hour — same cross-container staleness guard as
   `gb_download_cron.sh`, same caveat: **comment out the crontab line before
   running the pipeline by hand**, the guard is a backstop, not a
   substitute.
2. Compares `tar_read(gb_release, store = "_targets_gb_store")` against
   `tar_read(gb_release, store = "_targets")` — no new marker file, reuses
   targets both stores already carry. Equal → exit 0 (the common case).
3. Diverge → `bash run.sh` (reads `IMAGE_TAG` out of `run.sh` itself rather
   than duplicating it), waits for the detached container to exit by
   polling `docker ps`, then checks `tar_meta(fields = "error")` for a clean
   finish.
4. Clean finish → runs `R/publish_figshare.R` and `R/snapshot_ftol_data.R` in
   one container, as the host user (`HOST_UID`/`HOST_GID`, so nothing in the
   repo ends up root-owned) with the gitconfig mounted read-only at
   `/etc/gitconfig_persisted` (commit author only). The container has no `gh`
   and no credentials, so the snapshot **commits but does not push**; the
   script then pushes `ftol_data` from the host
   (`git -C ftol_data push origin main`), which has the `gh` credentials.
5. Either way, emails `joelnitta@gmail.com` via `send_release_ready_email()`
   — `"ready"` (snapshot + FigShare succeeded, your turn to review/release),
   `"publish_failed"` (pipeline succeeded but the snapshot/FigShare step
   errored), or `"pipeline_failed"` (`run.sh` itself errored).

## Enabling and disabling the job

Same recipe as `gb_download_cron.sh` (see `docs/nittalab_gb_download.md`),
substituting the schedule line:

```bash
# Disable
crontab -l | sed 's|^0 3 \* \* \*|#DISABLED# 0 3 * * *|' | crontab -
# Re-enable
crontab -l | sed 's|^#DISABLED# ||' | crontab -
```

## What to expect

- Most days: exits fast, "gb_release unchanged," nothing built — the
  `gb_download` project rarely has anything new for the `main` pipeline to
  pick up.
- When it does run: expect the same order-of-a-day-or-two runtime as a
  manual `run.sh` invocation (the `main` pipeline itself doesn't get faster
  by being cron-triggered), then a few more minutes for the publish step.
- An email always goes out on the "new data found" path — success, or
  either failure mode — so a silent multi-day gap here means it's still on
  the "unchanged" fast path, not stuck.
