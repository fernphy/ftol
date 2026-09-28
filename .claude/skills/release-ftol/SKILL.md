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
  current `ftol` HEAD, snapshotting is done.
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
   deposit `19474316` via `upload_to_figshare()`; `overwrite = TRUE`, no
   manual delete step).
8. Commit any final `ftol` code changes yourself, as usual.
9. `Rscript R/snapshot_ftol_data.R` — see `run-tar-make`'s "After a clean
   run" section for its five checks. Safe to run non-interactively once
   they're green.

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

## Out of scope

- Step 17 (Mastodon) — not covered by this skill; do it yourself.
- First-digit version bumps anywhere — always ask, never infer.
