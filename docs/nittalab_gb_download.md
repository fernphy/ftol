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
   **TODO for nittalab setup: confirm the actual mount path for the external
   drive Joel has been manually archiving to, and use that path below.** If
   this isn't findable/known, ask before guessing — silently pointing this
   at the wrong path would mean archival quietly does nothing.

## Host crontab entry

There's roughly a 3-month gap between GenBank releases, so a low-frequency
check is enough — daily is a safe, simple default. Since a full run can take
1-2 days, guard against overlapping invocations with `flock`; `targets`
itself also refuses to run twice against the same locked store, so a naive
overlap is harmless (fails fast) but `flock` keeps the logs clean:

```cron
0 0 * * * flock -n /path/to/ftol/.gb_download.lock -c '\
  cd /path/to/ftol && \
  docker run --rm \
    -v $(pwd):/wd -w /wd \
    -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) \
    -e GB_DL_ARCHIVE_DIR=/path/to/external/drive \
    joelnitta/ftol:latest \
    Rscript -e "Sys.setenv(TAR_PROJECT = \"gb_download\"); targets::tar_make(script = \"_targets_gb.R\")" \
  >> /path/to/ftol/logs/gb_download_cron.log 2>&1'
```

Adjust the image tag to match whatever `run.sh` currently uses. Create
`logs/` first if it doesn't exist (it's gitignored).

## What to expect

- Most days, this will do nothing: `_targets_gb.R`'s `release_check` target
  errors out immediately with "No new GenBank data available" whenever the
  latest NCBI release matches what's already installed. That's normal, not a
  failure to act on.
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
