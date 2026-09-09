# Setting up automated GenBank downloads on the nittalab server

This is a guide for Claude (or a human) implementing recurring, automated
GenBank release checks on the nittalab server, following a migration away
from the old `R/setup_gb.R` script (deleted) to a `targets` pipeline,
`_targets_gb.R`.

## Why this changed

The old script downloaded the entire GenBank plant division (~1.24 TB
compressed) to `scratch/` before filtering it down to ferns. `_targets_gb.R`
streams it instead: each of the ~3,343 plant-division files is downloaded,
scanned with a fast pre-check, parsed and filtered to target accessions, and
deleted individually, so peak disk use is ~1 file at a time. This was
originally driven by a devcontainer with far less disk than the division
needs, but the resulting pipeline is a straightforward improvement
everywhere: dynamic branching over per-file targets gives free resumability
(a crashed/interrupted run just picks up where it left off on the next
`tar_make()` invocation, no custom progress-tracking needed) and safe modest
parallelism (each per-file branch is a pure function returning data, not a
shared-file mutation, so concurrent workers don't conflict).

Cron scheduling also moved: it used to be baked into the Docker image
(`Dockerfile`'s now-removed "Cron" section, running `cron -f` inside a
long-lived container). That's gone. Baking "when to run" into the image
conflated it with "what the image does," which breaks down for anything that
isn't a long-lived container (devcontainers, CI, image rebuilds). Scheduling
now lives at the host level, calling `docker run` on demand — the same
pattern `run.sh` already uses for the main FTOL pipeline.

## One-time setup

1. Confirm the FTOL repo is checked out at a stable path on the host (not
   inside a transient container), e.g. `/path/to/ftol`.
2. Confirm `.secrets/` (gmailr OAuth credentials) is present in the repo root
   — needed for the start/finish notification emails to
   `joelnitta@gmail.com`. These are already gitignored; copy them over
   separately if setting this up on a fresh checkout.
3. Confirm the `joelnitta/ftol` Docker image is built and available (see
   `run.sh` for the current tag).
4. **Fill in the external archive path.** `_targets_gb.R` supports an
   optional `GB_DL_ARCHIVE_DIR` environment variable: if set to a directory
   that exists, the pipeline automatically copies the outgoing release's
   database (before overwriting it) into `<archive_dir>/gb_release_<N>/`,
   matching the existing manual `gb_release_<N>/` naming convention. This
   matters because NCBI's FTP server only serves the *current* release's
   flatfiles — once a release is superseded, a filtered database built from
   the old one can never be regenerated. A single local `.bak` copy is
   always kept regardless (zero configuration needed for that much), but the
   external drive is where multiple past releases should live long-term.
   **nittalab setup: the archive location is
   `/mnt/jnitta/project_data/ftol_genbank_raw`** (NFS mount from the Synology
   NAS), which already holds the manually-created `gb_release_<N>/` dirs. The
   pipeline writes new `gb_release_<N>/` subdirs there. Because the pipeline
   runs inside the container, this host path must be bind-mounted in and
   `GB_DL_ARCHIVE_DIR` set to the in-container path — the crontab entry below
   does both (`-v /mnt/jnitta/project_data/ftol_genbank_raw:/archive`,
   `-e GB_DL_ARCHIVE_DIR=/archive`).

## Host crontab entry

There's roughly a 3-month gap between GenBank releases, so a low-frequency
check is enough — daily is a safe, simple default. Since a full run can take
1-2 days, overlapping invocations must be prevented:

- **`flock -n` in the crontab** stops one cron tick from starting while the
  previous one is still running.
- **A staleness check in the wrapper** stops a tick from colliding with a
  run started *outside* cron — e.g. a manual `tar_make` in a devcontainer.
  This matters because `targets`' own "a pipeline is already running" guard
  is **PID-based and does not work across containers** (separate PID
  namespaces — the cron container cannot see the devcontainer's process, and
  vice versa). Two `tar_make` processes on one store is not automatically
  refused. The wrapper checks the mtime of `_targets_gb_store/meta/progress`
  (the running pipeline touches it constantly) and skips the tick, exiting 0,
  if it was modified in the last hour.

Still: **before running the pipeline by hand, comment out the crontab line**,
and re-enable it afterward. The staleness check is a backstop, not a
substitute.

The `docker run` invocation lives in a small wrapper script,
[`gb_download_cron.sh`](../gb_download_cron.sh) (repo root, committed
alongside `run.sh`), so the crontab line stays readable and the quoting
stays sane. It brackets the run with `=== gb_download cron start/end ===`
markers, applies the staleness check above, runs `targets::tar_make()` in the
container, and on any non-zero exit prints the recorded target errors
(`tar_meta(fields = "error")`) so the log says *why* it stopped.

Key pieces of that `docker run`:

- `-v /home/jnitta/ftol:/wd -w /wd` — the repo is the working directory, same
  as `run.sh`.
- `-e HOST_UID=$(id -u) -e HOST_GID=$(id -g)` — `entrypoint.sh` then runs R
  as `jnitta` rather than root, so files written into the bind mount (and the
  archive) are owned correctly.
- `-e TAR_PROJECT=gb_download` — `targets` picks up `script: _targets_gb.R`,
  `store: _targets_gb_store` and `reporter_make: summary` from
  `_targets.yaml`. The `summary` reporter keeps the log to a single
  self-rewriting status line; a full run maps over ~3,343 per-file branches,
  so `verbose` (the default) would write thousands of lines.
- `-v /mnt/jnitta/project_data/ftol_genbank_raw:/archive` +
  `-e GB_DL_ARCHIVE_DIR=/archive` — the outgoing release archive (see above).

jnitta's crontab entry (daily at 00:00):

```cron
0 0 * * * /usr/bin/flock -n /home/jnitta/ftol/.gb_download.lock /home/jnitta/ftol/gb_download_cron.sh >> /home/jnitta/ftol/logs/gb_download_cron.log 2>&1
```

`logs/` must exist and be writable by `jnitta` (it's gitignored). `flock`
creates `.gb_download.lock` on first run (also gitignored). Adjust the image
tag if `run.sh` moves off `joelnitta/ftol:latest`.

## Enabling and disabling the job

`crontab -e` is interactive; to toggle the entry from a script or a
non-interactive session, filter the crontab through `crontab -`.

Disable (comment the schedule line out, keeping it in place):

```bash
crontab -l | sed 's|^0 0 \* \* \*|#DISABLED# 0 0 * * *|' | crontab -
```

Re-enable:

```bash
crontab -l | sed 's|^#DISABLED# ||' | crontab -
```

Confirm with `crontab -l`. Disable the job before running the pipeline by
hand (see the overlap note above), and whenever the server is being worked
on; re-enable it afterward.

## What to expect

- Most days, this will do nothing: `_targets_gb.R`'s `release_check` target
  errors out immediately whenever the latest NCBI release matches what's
  already installed, so the run exits non-zero with `errored | 1` in the
  summary line and nothing built. That's normal, not a failure to act on.
  Confirm the reason with
  `targets::tar_meta(fields = "error", complete_only = TRUE)` (expect
  "No new GenBank data available").
- When a new release does appear, expect the full run to take on the order
  of a day or two (network-bound: NCBI serves the whole plant division
  either way, ~1.24 TB compressed, whether streamed file-by-file or
  downloaded in bulk). An email goes out to `joelnitta@gmail.com` when it
  starts and when it finishes.
- If a run gets interrupted (host reboot, network drop, container killed
  mid-way), the next cron tick just resumes — completed per-file branches
  are cached and skipped, only unfinished ones re-run.
- The official `_targets/user/data_raw/restez/gb_release.txt` (the file
  `release_check` reads to decide whether there's new data) is deliberately
  the *very last* thing the pipeline touches, after the database, README,
  FigShare archive, and taxdmp are all confirmed built — so a crash at any
  point beforehand leaves that gate file untouched and a resumed run
  correctly picks the remaining work back up, rather than concluding
  (wrongly) that nothing is left to do.
- After a successful run: follow the FigShare/release-version steps already
  documented in `docs/updating.md` from there.

## Debugging

- `_targets_gb.R` accepts a handful of environment variable overrides
  originally added for dry-run validation, all safe to leave unset in
  production (they default to normal behavior): `GB_DL_DATA_RAW`,
  `GB_DL_SCRATCH`, `GB_DL_SEND_EMAIL`, `GB_DL_FILE_CAP` (caps the number of
  plant-division files processed — useful for a quick smoke test against a
  handful of real files before trusting a fresh setup with the real
  multi-day run), `GB_DL_ARCHIVE_DIR` (above).
- Inspect progress/failures with `targets::tar_manifest()`,
  `targets::tar_meta(fields = error, complete_only = TRUE)`, or
  `targets::tar_progress()` against the `_targets_gb_store` store (set
  `Sys.setenv(TAR_PROJECT = "gb_download")` first, or pass
  `store = "_targets_gb_store"` directly).
