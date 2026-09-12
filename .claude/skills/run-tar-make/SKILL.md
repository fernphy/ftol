---
name: run-tar-make
description: Launch `targets::tar_make()` for the FTOL pipeline — a full run or up to one specific target — log it consistently under `logs/`, and monitor it through to completion or the next error. Use whenever asked to run/kick off/relaunch/resume the pipeline, run tar_make, "run it up to `<target>`", or to check on a pipeline that's already running. Covers the log naming convention (shared with `run.sh`), launching in the background, monitoring progress, and recovering from a stale `_targets/meta/process` lock or a leftover zombie process. For diagnosing *why* a run fails on a taxonomy-related target, use the `ppg-taxonomy-update` skill instead — this skill is purely about the mechanics of launching and watching a run.
user-invocable: true
---

# run-tar-make — launch and monitor the FTOL pipeline

## Two execution contexts — figure out which one you're in first

This pipeline gets run from two genuinely different places, and they call for different
behavior:

- **Interactive dev container** — this environment (or an equivalent devcontainer session
  on nittalab), with Claude or Joel driving `tar_make()` directly for development: testing
  a fix, checking one target, iterating on a taxonomy update, etc. Detect it with
  `[ -f /.dockerenv ] && [ -z "$IMAGE_TAG" ]` — inside *a* container, but not one `run.sh`
  launched for this run. In this context: just call `tar_make()` (or
  `tar_make(names = "<target>")`) directly as below, in the background, logged under
  `logs/` — no need to touch Docker at all, you're already inside the image. This is
  **never** the pathway for producing an official release: `image_tag` records a
  placeholder here (see below), and `snapshot_ftol_data.R` refuses to publish a snapshot
  built from one.
- **nittalab server, bare host, via `run.sh`** — the production pathway. `run.sh` launches
  a *fresh* container specifically for that one pipeline run and passes it `IMAGE_TAG`, so
  the resulting `image_tag` target is trustworthy provenance for a release (see
  `docs/updating.md`). Detect it by the *absence* of `/.dockerenv` — a bare shell on the
  host, not already inside any container. In this context, if asked to run the real/
  official/production pipeline (or to prep for a data release), invoke `bash run.sh` and
  read its logging output — do **not** hand-roll an equivalent `docker run
  ... Rscript -e 'targets::tar_make()'` call; `run.sh` already sets `IMAGE_TAG` and the
  `logs/` convention correctly, and reimplementing it risks getting one of those wrong.

The rest of this skill (logging, launching, monitoring) applies the same way in the dev
container case; for the `run.sh` case, `run.sh` already handles logging/launching itself —
this skill's job there is mainly monitoring the log it produces and recovering from a
stale lock or zombie if something got interrupted.

## Logs: always under `logs/`, one naming scheme

**Every `tar_make()` log lives in `logs/`, named `logs/tar_make_<YYYYMMDD_HHMMSS>.log`.**
`logs/tar_make_latest.log` always symlinks to the most recent one, so "what's the current
log" never requires globbing dates. `run.sh` (the container/Docker entry point) follows the
same convention — this is the one place to look regardless of how the run was started.
`logs/` also holds unrelated `gb_download_*`/`gb_fix_*` logs from the separate GenBank
download pipeline (`_targets_gb.R`) — leave those alone, they're a different naming scheme
for a different pipeline.

Before starting a run:

```bash
mkdir -p logs
LOG="logs/tar_make_$(date +%Y%m%d_%H%M%S).log"
ln -sf "$(basename "$LOG")" logs/tar_make_latest.log
```

## Launching

Run from the project root (`/home/rstudio/ftol/ftol`) so `renv` is picked up, in the
background so it survives past this turn:

```bash
nohup Rscript -e 'targets::tar_make(reporter = "verbose_positives")' > "$LOG" 2>&1 &
```

No target needs excluding for a dev-container run: `image_tag` (docker image provenance,
used only when publishing a release — see `get_docker_tag()`) returns a placeholder rather
than erroring when `IMAGE_TAG` isn't set, so a plain `tar_make()` always succeeds here. That
env var is only ever set by `run.sh`, which launches a fresh container specifically for
that run.

For a **narrower run** — e.g. rebuilding just far enough to check one target after a fix —
pass `names = "<target>"`; `tar_make()` still builds every upstream dependency, it just
stops once that target is done:

```bash
nohup Rscript -e 'targets::tar_make(names = "non_mono_check",
  reporter = "verbose_positives")' > "$LOG" 2>&1 &
```

## Monitoring

Name resolution / taxonomy targets finish within minutes; alignment and tree-building are
the long compute phase (hours+). `tar_make()` stops dispatching new targets at the first
error but lets already-running ones finish, so the `_targets` store stays safe to inspect
mid-run (`tar_meta(fields = "error")`, `tar_progress()`, `tar_workspace(<target>)`).

- Quick look: `tail -f logs/tar_make_latest.log`, or grep it for
  `completed target|errored|Execution halted`.
- For a standing watch across turns, use the `Monitor` tool on
  `tail -n 0 -f logs/tar_make_latest.log | grep -E --line-buffered "..."` — filter tightly
  (errors, `Execution halted`, and only the milestone targets you care about) since the
  full target list is chatty; a filter that's too broad gets auto-throttled.
- To know when the run itself ends (success or error), a background `Bash` loop works
  better than watching the log, since state can go quiet without the process exiting:
  ```bash
  while pids=$(pgrep -f 'targets::tar_make'); do
    alive=0
    for p in $pids; do
      st=$(ps -o stat= -p "$p" 2>/dev/null | tr -d ' ')
      case "$st" in Z|"") ;; *) alive=1;; esac
    done
    [ "$alive" -eq 0 ] && break
    sleep 30
  done
  echo "tar_make finished"
  ```
  Match on `stat` rather than plain existence — see the zombie note below.

## Recovering from an interrupted run

- **Stale process lock**: `tar_make()` refuses to start ("Process ID … is already running
  a {targets} pipeline") if `_targets/meta/process` names a pid from a run that got killed
  or whose container/session ended without cleanup. Check whether that pid is actually
  alive (`ps -p <pid>`); if not, `rm -f _targets/meta/process` before relaunching.
- **Zombie processes**: this container's PID 1 is `sleep infinity`, which does not reap
  children, so a killed run's process can sit as a `Z` (defunct) zombie indefinitely — it's
  harmless (holds no real resources) but `tail --pid=<pid>` and `kill -0 <pid>` both still
  "succeed" against a zombie, so waiters built on those hang forever. Match on
  `pgrep -f 'targets::tar_make'` plus `ps -o stat=` not equal to `Z` instead (as in the
  monitoring loop above).
- To actually stop a run: kill the outer `Rscript` pid and its `crew`/`mirai` children
  (`pkill -f 'mirai::dispatcher'`, `pkill -f 'crew::crew_worker'`), then clear the process
  lock as above before relaunching.

## Fixed hang to know about: `mpcheck_monophy`

Historically `check_monophy()` (target `mpcheck_monophy`) hung indefinitely — one core
pegged, no progress for hours, "fixed" only by re-running. Cause: it ran with `deployment =
"main"` (executing in the `tar_make` controller process itself) and spun up a
`future::multisession` PSOCK cluster there, which deadlocks against `crew`'s `mirai` event
loop sharing that process; `deployment = "main"` also froze every other crew branch
meanwhile. Fixed 2026-09 by removing `deployment = "main"` and the internal
`future`/`furrr` (now plain `purrr::map_lgl` — `ape::is.monophyletic()` is ~0.02–0.2 s/call,
so a whole locus is seconds to ~3 min, and `crew` parallelizes across the 7 loci). If this
exact hang reappears on an old checkout, that's the fix to reapply. More generally: **never
nest `future::multisession` inside a `deployment = "main"` target** — let `crew` do the
parallelism instead. The other `future::multisession` calls in the codebase
(`parse_gb_gene` / `parse_gb_spacer` / `extract_ncbi_names`) run inside isolated crew
workers rather than the controller, so they're lower-risk, but are candidates for the same
fix if one of them ever hangs too.
