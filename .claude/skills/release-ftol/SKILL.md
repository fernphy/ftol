---
name: release-ftol
description: Walk through publishing an FTOL update across all 5 repos (docs/updating.md steps 1-16 -- step 17, Mastodon, is out of scope), running the automated prep at each step and pausing for explicit confirmation before anything public/irreversible (a GitHub release, rsconnect::deployApp(), or a push to fernphy.github.io). Use when asked to release/publish/ship an FTOL update, run through docs/updating.md, or continue a release after main_pipeline_cron.sh's "ready" email.
user-invocable: true
---

# release-ftol — publish an FTOL update end to end

Companion to `run-tar-make` (pipeline mechanics) and `ppg-taxonomy-update`
(taxonomy failures) — this skill is the release orchestration layer on top.
Read `docs/updating.md` for the full narrative; this skill is the "what to
actually run, and where the confirmation gates are" version.

## Where things stand when you start

Check which of these already happened before doing anything:

- Has `tar_make()` finished cleanly? `tar_meta(fields = "error", complete_only
  = TRUE)` should have zero rows. If not, this isn't ready yet — that's
  `run-tar-make`'s job, not this skill's.
- Has `R/publish_figshare.R` / `R/snapshot_ftol_data.R` already run? Check
  `git -C ftol_data log -1` — if its message has `code=<hash>` matching the
  current `ftol` HEAD, snapshotting is done. For FigShare, don't trust a
  log or handoff note: compare the *public* record (unauthenticated API)
  against the local files; "uploaded and verified" only means the pending
  edit is staged (see steps 5-7).
- If `main_pipeline_cron.sh` (see `docs/nittalab_main_pipeline_cron.md`) is
  set up on the host, both of the above may already be done and you're
  picking this up from its "ready for version bump/push" email — in that
  case, start directly at step 10 below.
- If Claude starts a session and finds the pipeline finished and waiting at
  a gate with no indication the user has seen it, use `PushNotification` to
  say so — best-effort, not the primary trigger (that's the cron chain).

## Steps 1-9 (in-repo, already automated)

1. GenBank download — cron'd (`gb_download_cron.sh`,
   `docs/nittalab_gb_download.md`).
2. Age sanity check — now a target, `age_comparison` (`results/
   age_comparison.csv` under `_targets/user/results/`). Skim it for large
   `diff` values before proceeding; not a hard gate, ages can legitimately
   shift.
3. Change log — now automatic (`update_changelog()`, `R/functions.R`):
   `reports/{input_data_readme,ftol_data_readme}/changelog_history.txt` get
   a new dated entry when `gb_release` advances. Only touch these by hand
   for a non-routine entry (e.g. "Add taxdmp.zip").
4. `bash run.sh` — see `run-tar-make`.
5-7. FigShare publish — `Rscript R/publish_figshare.R` (uploads
   `restez_sql_db.tar.gz`, `taxdmp.zip`, `README.genbank`, `README.txt` to
   deposit `19474316` via `upload_to_figshare_verified()`; overwrites, no
   manual delete step). **A printed "verified by checksum" does NOT mean it
   is public.** On an already-public deposit, uploads only change the
   owner's *pending edit* (https://figshare.com/account/articles/19474316);
   the public record keeps serving the old files until the owner clicks
   Publish. Procedure:

   1. Run it in tmux/detached (the restez upload is large).
   2. List the pending files:
      `curl -s -H "Authorization: token $FIGSHARE_TOKEN"
      https://api.figshare.com/v2/account/articles/19474316/files`
      (token is in `.Renviron`; never echo it). Retried uploads leave
      same-named duplicates with identical md5s. Remove them with
      `bash figshare_delete_files.sh <id>:<name>:<md5> ...` — it only
      deletes an id if its name and md5 still match, never publishes, and
      is pre-approved in `.claude/settings.local.json`. Keep one of each
      name (`ref_aln.tar.gz` is untouched and must remain).
   3. → **Confirmation gate: the user publishes.** Tell them the pending
      files look right (names, sizes, md5s equal to the local files under
      `_targets/user/data_raw/`) and ask them to click Publish. Do not call
      the publish API yourself.
   4. Verify the *public* record afterwards: unauthenticated
      `GET https://api.figshare.com/v2/articles/19474316` shows a bumped
      `version`, and the downloaded `README.genbank` names the new GenBank
      release. Only then continue.
8. Commit any final `ftol` code changes yourself, as usual.
9. Snapshot — do the "Before the snapshot" preflight below, then
   `Rscript R/snapshot_ftol_data.R` inside the container (see
   `run-tar-make`'s "After a clean run" section for the command and its
   checks). Safe to run non-interactively once the preflight is green.

### Before the snapshot (preflight; each item bit us once)

- **`restez_sql_db_hash` / `ref_aln_hash` in `R/snapshot_ftol_data.R`**
  are hard-coded and go stale whenever the archives are rebuilt. Don't just
  paste a new hash: first confirm the local archive is what FigShare serves
  (the `gb_release.txt` inside `restez_sql_db.tar.gz` and `README.genbank`
  should state the release you're shipping, and the md5s should match the
  public record). Then update the hash (`contentid::content_id(path)`) and
  have the user commit/push — the script requires a clean code repo.
- **Code repo must be clean as seen from the container**, which doesn't
  read the host's global git ignore. Untracked local-only files (e.g.
  `.claude/settings.local.json`) make the check fail; list them in the
  repo's `.git/info/exclude`, not just the global ignore.
- **`ftol_data` must be able to fast-forward.** `git -C ftol_data status -sb`
  — if it is behind origin and has local edits to tracked non-data files
  (typically `release.R`, left modified by an earlier session), `git_pull`
  fails with "1 conflict prevents checkout". Fix: `git -C ftol_data stash
  push -- release.R` before the snapshot (do NOT `git checkout --` it; that
  discards work and is blocked), then after the snapshot `git stash pop`,
  keep the automated version, and have the user commit/push `release.R`
  (step 10 also needs a clean `ftol_data`).
- If a root-owned file still shows up (`find ~/ftol -user root`), it came
  from a container run without `HOST_UID`/`HOST_GID` or from a dev
  container session started before `remoteUser` was set; fix once with
  `sudo chown -R jnitta:jnitta ~/ftol`.
- **Run the snapshot/publish container as the host user, and push from the
  host.** Use the `docker run` in `run-tar-make`'s "After a clean run"
  section: `-e HOST_UID=$(id -u) -e HOST_GID=$(id -g)` (no root-owned
  files left in `ftol_data/.git`, issue #38) and the gitconfig mounted at
  `/etc/gitconfig_persisted`. The image has no `gh` and the container has
  no credentials, so the snapshot commits and prints a reminder instead of
  pushing -- then run `git -C ftol_data push origin main` on the host.
  Docker is only for the image-bound scripts (`tar_make()`, FigShare
  publish, snapshot); everything else runs on the host.
- `write_cc0()` doesn't need the contentid cache (a `docker run --rm`
  starts with an empty one) as long as `ftol_data/LICENSE` is already the
  CC0 text.

## Steps 10-16 (cross-repo, confirmation-gated)

Cloned as siblings of `ftol_data` (nested inside this repo, gitignored):
`ftolr`, `ftol_vis`, `ftol_shiny`, `fernphy.github.io`. **Strict order** —
each step's script depends on the previous one's release actually existing
(a download-by-tag that 404s is the ordering guard, not a stop-you-first
check):

**10. `ftol_data` release.** `cd ftol_data && Rscript release.R` — derives
`new_ver`/`notes` from the `ftol` pipeline's `gb_release`/`date_cutoff`
targets and the last release's notes (semver rule: 3rd digit if GenBank
release is unchanged, 2nd digit if new; **1st digit is never inferred** —
flag it and ask if you think a breaking change is warranted). Commits/pushes
the CFF bump automatically, then prints the `gh release create` command
instead of running it.

→ **Confirmation gate.** Show the user the proposed `new_ver`/`notes` and
the exact command. Only run it after they say go.

**11-12. `ftolr` data update.** `cd ftolr && Rscript update_data_ver.R`
— rewrites `ft_data_ver.R`'s three literal values, runs
`data-raw/import_data.R`, commits. Fully mechanical, safe to run once step
10's release exists.

**13. `ftolr` release.** `Rscript inst/release.R` — already derives its own
`new_ver`/notes from `ft_data_ver()` (untouched, already correct). Its last
two lines push and call `gh release create`.

→ **Confirmation gate** before running this script (or at least before its
final `gh release create` line) — it's a public release affecting every
downstream `ftolr` install.

**14. `ftol_vis`.** `cd ftol_vis && Rscript update_ftolr.R` — `renv::install`
+ `tar_make()`. Purely local, safe to run unattended.

**15. `ftol_shiny`.** `cd ftol_shiny && Rscript update_ftolr.R` —
`renv::install`, snapshot, commit, smoke-test via `shiny::runApp()`. Purely
local up through the smoke test.

→ **Confirmation gate** before `rsconnect::deployApp('ftol_explorer')` (the
script prints this command rather than running it) — goes live publicly.
Needs `rsconnect` installed in `ftol_shiny`'s renv first if it isn't
already.

**16. `fernphy.github.io`.** `cd fernphy.github.io && Rscript update_ftolr.R`
— freezes the outgoing "Current" version into `downloads.Rmd`'s "Past
versions" (needs the *old* ftolr still installed, hence this runs before
step 14's `renv::install` equivalent here), drafts a `news.Qmd` entry. Read
over the drafted headline wording — that part is a reasonable default, not
polished copy. Then `renv::install("fernphy/ftolr")` and
`quarto render` (not `quarto preview` — that only rebuilds pages whose
*source* changed, and misses content that only changed because `ftolr`
changed). Needs `quarto` installed in this container, or run this step from
wherever it is.

→ **Confirmation gate** before pushing — GH Actions deploys the site
immediately on push to `main`.

## Running steps 10-16 on the bare host (not the dev container)

The scripts' hard-coded sibling paths (`../ftol`) assume the dev-container
layout; on the host the repos are nested (`~/ftol/{ftol_data,ftolr,...}`),
so the pipeline's `_targets` is at `..`. Each bit us once:

- `ftol_data/release.R` and `ftolr/update_data_ver.R` need
  `ftol_repo <- if (dir.exists("../ftol/_targets")) "../ftol" else ".."`
  (done in `release.R`, committed; `update_data_ver.R` is untracked in
  `ftolr`, patched in place). `release.R` also needs a dashed cutoff date
  (`gsub("/", "-", ...)`), as does `update_data_ver.R`, otherwise
  `ft_data_ver("cutoff")` becomes `2026/08/02` and ends up in release notes.
- The host `gh` is old (2.4.0): `gh release create ... --latest` fails with
  `unknown flag`. Omit `--latest` (the newest non-prerelease becomes Latest
  anyway). `ftolr/inst/release.R` already omits it.
- `ftolr`: untracked `update_data_ver.R` makes `inst/release.R`'s
  clean-repo check fail and would be swept into the script's own commit --
  add it to `ftolr/.git/info/exclude` and exclude it from the commit.
  `devtools::build_readme()` in `inst/release.R` always changes `README.md`
  (example output) -- commit it ("Update README") before the push/release.
  `git fetch --tags` in `ftol_data` first, or `git describe` misses the new
  tag.
- Sibling repos (`ftol_vis`, `ftol_shiny`, `fernphy.github.io`) have empty
  renv libraries on the host: run `renv::restore(prompt = FALSE)` in tmux
  first (~15 min, compiles from source; most is cached afterwards). Needs
  system libs `libharfbuzz-dev libfribidi-dev libtiff-dev libjpeg-dev
  libwebp-dev`. `ftol_shiny`'s `update_ftolr.R` ends with a blocking
  `shiny::runApp()` and needs `gert` (not in its renv) -- run restore +
  `renv::install("fernphy/ftolr")` + `renv::snapshot()` by hand, commit
  `renv.lock`, and smoke-test headless (`runApp(..., launch.browser =
  FALSE)` + `curl`).
- `fernphy.github.io` renv is 0.15.5 with `renv.config.pak.enabled = TRUE`
  in `.Rprofile`; `renv::install()` then fails with "Cannot parse package".
  Run it with `options(renv.config.pak.enabled = FALSE)`. The restore must
  happen BEFORE `update_ftolr.R` (it needs the old `ftolr` installed).
  `update_ftolr.R` doesn't update the previous release's news link: change
  `#current-v<old>` to `#v<old>` in `news.Qmd` and delete its FIXME.
- `renv::snapshot()` on the host rewrites many `Repository` fields
  (`CRAN` -> `RSPM`/PPM URL) and reorders entries. Rebuild `renv.lock` from
  `git show HEAD:renv.lock` with only the `ftolr` Version/RemoteSha/Hash
  lines changed, so the commit is a 3-line diff.
- `git add` refuses paths in `.gitignore` even when tracked (`ftol_vis`'s
  `_targets/user/taxonium/*`): use `git add -u <path>`.
- shinyapps.io deploy needs an rsconnect account on the host
  (`rsconnect::setAccountInfo()`, done once by the user in `ftol_shiny`);
  target the existing app explicitly: `rsconnect::deployApp("ftol_explorer",
  appName = "ftol_explorer", account = "fernphy", server = "shinyapps.io",
  forceUpdate = TRUE)`. The auto-mode classifier blocks this and the
  `gh release create` calls even with an allow rule -- the user runs them,
  or approves the prompt.

## Working agreements that applied to this release

- **Commits only on request** (org policy in `/etc/claude-code/CLAUDE.md`),
  including the scripts' own commits (`release.R`, `update_data_ver.R`,
  the `ftolr` release prep). Get one explicit OK for "let the script make
  its N local commits" per step, and otherwise ask the user to commit/push
  code changes (the snapshot requires a clean, pushed code repo).
- **Public actions are the user's call even when approved in chat.** The
  auto-mode classifier blocked `gh release create` (once), the shiny
  deploy, and FigShare API deletes, regardless of allow rules for the
  exact command. Don't retry or rephrase: hand the user the exact command.
  Narrow, scripted helpers (`figshare_delete_files.sh`) were allowed.
- **Verify the public result after every public step**, don't trust the
  script's last line: `gh release list -R fernphy/<repo> -L 1`; for the
  site `gh run list -R fernphy/fernphy.github.io -L 1` then `curl` the page
  for "Current: v<new>"; for the app `curl` for HTTP 200 and the version
  string; for FigShare the unauthenticated API (see steps 5-7).
- **The host does have R 4.5 + renv** (the "no local R on the bare host"
  note in `run-tar-make` is about the *pipeline*, which needs the image).
  Steps 10-16 run on the host; only the image-bound scripts (`tar_make()`,
  FigShare publish, snapshot) run in Docker, as the host user.
- Run anything slow (renv restores, `renv::install`, FigShare upload) in
  detached tmux with a log under `logs/` (user's global CLAUDE.md), and
  put R snippets in a script file rather than nesting quotes in
  `tmux new-session "..."` (a mangled `Rscript -e` caused a spurious
  "Cannot parse package" error).
- Delete `RELEASE_HANDOFF.md` and drop leftover `git stash` entries when
  the release is finished; remove one-off permission rules from
  `.claude/settings.local.json`.

## Out of scope

- Step 17 (Mastodon) — not covered by this skill; do it yourself.
- First-digit version bumps anywhere — always ask, never infer.
