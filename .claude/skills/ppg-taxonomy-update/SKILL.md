---
name: ppg-taxonomy-update
description: Fix the taxonomy / name-resolution errors that block `targets::tar_make()` after the PPG taxonomic database is bumped to a new version in `_targets.R` (`load_ppg(ver = ...)`). Covers stale workarounds in `modify_ppg()`, broken `match` targets in `pterido_manual_match.csv`, `unchecked` names needing `format_ppg_for_ts()` `recs_keep` exceptions, the inline `plastome_manual_match` list, the `create_patel_inclusion_list()` name repairs, and the `pterido_names_to_inspect` / `ncbi_accepted_names_map` human-review checkpoint. Use when a GenBank/PPG update makes the pipeline error in `modify_ppg`, `ppg_db`, `ppg_ref_names`, `match_results_raw_round_*`, `plastome_metadata_renamed`, `patel_inclusion_list`, or `ncbi_accepted_names_map`.
user-invocable: true
---

# ppg-taxonomy-update — clear taxonomy errors after a PPG version bump

Each new GenBank release goes with a newer PPG (Pteridophyte Phylogeny Group) taxonomy.
Bumping `ppg_raw = load_ppg(ver = "0.0.0.9xxx")` in `_targets.R` rebuilds `ppg_full`,
`ppg_db`, `ppg_ref_names` and every name-resolution target, and that reliably surfaces the
same handful of failure modes: workarounds baked into the code for problems PPG has since
fixed, and `pterido_manual_match.csv` rows whose target name PPG has since renamed or
removed. This skill walks each one in the order the pipeline hits it. `docs/updating.md`
is the surrounding release checklist; this is only the "make `tar_make()` get past name
resolution" part. For actually launching/monitoring a `tar_make()` run (log location,
background launch, recovering from a stale lock or zombie process) use the `run-tar-make`
skill — it's a separate concern from diagnosing *why* a run failed.

## Orientation — the data flow

```
load_ppg(ver)            download ppg.csv for that tag  → ppg_raw
  │
modify_ppg()             hand-patch ppg_raw: drop duplicate taxonIDs, dct_add_row /
  │                      dct_modify_row for taxa missing/wrong in PPG        → ppg_full
  │
format_ppg_for_ts()      drop names that aren't accepted/synonym or have a bad
  │                      nomenclaturalStatus; `recs_keep` is the taxonID allow-list
  │                      that survives that filter anyway                    → ppg_db
  │
ts_parse_names(unique(ppg_db$scientificName))                                → ppg_ref_names
  │
ts_match_names(query = <NCBI names>, reference = ppg_ref_names,
  │            manual_match = manual_matches)   ← manual_matches = pterido_manual_match.csv
  │            → match_results_raw_round_1/2/3 → …_resolved_… → match_results_resolved_all
  │
inspect_ts_results()     collect every auto_fuzzy / no_match name            → pterido_names_to_inspect
  │
make_ncbi_accepted_names_map(strict = TRUE)   errors unless the above is 0 rows
```

Two other consumers of the same reference + CSV:

- `resolve_pterido_plastome_names()` (target `plastome_metadata_renamed`) has its **own**
  inline `plastome_manual_match` tibble that is `bind_rows()`'d with `pterido_manual_match.csv`.
- `create_patel_inclusion_list()` (target `patel_inclusion_list`) has a `case_when` block
  that rewrites species names from the 2019 Patel supplement to current PPG spelling before
  `ts_resolve_names()`.

## Get the new ppg.csv for inspection

A clone of the PPG repo lives at `/home/rstudio/ftol/ppg` (all release tags). Prefer it:

```bash
git -C /home/rstudio/ftol/ppg fetch --tags -q
git -C /home/rstudio/ftol/ppg show v0.0.0.9xxx:data/ppg.csv > /tmp/ppg_new.csv
```

Fallback if that clone is missing: `curl -sL -o /tmp/p.zip https://github.com/pteridogroup/ppg/archive/refs/tags/v0.0.0.9xxx.zip && unzip -p /tmp/p.zip 'ppg-0.0.0.9xxx/data/ppg.csv' > /tmp/ppg_new.csv`

`ppg.csv` columns: `taxonID, scientificName, scientificNameAuthorship, taxonRank,
parentNameUsageID, nomenclaturalStatus, namePublishedIn, taxonomicStatus,
acceptedNameUsageID, created, modified`. `taxonID` is the stable WFO id.

## Working method — fix standalone, then run the pipeline

Don't iterate by re-running `tar_make()` (minutes per loop, and it stops the whole run).
Reproduce each failing target's core logic in a throwaway `Rscript` that reads upstream
inputs with `tar_read()` and the freshly-built `ppg_db` / edited CSV, fix until it's clean,
then do one real pipeline run (see the `run-tar-make` skill for how to launch and log it).
Skeleton:

```r
suppressPackageStartupMessages({
  library(targets); library(tidyverse); library(taxastand); library(assertr)
  library(readr); library(tidyr); library(dwctaxon)
})
source("R/functions.R")
ppg_raw  <- load_ppg(ver = "0.0.0.9xxx")     # or readr::read_csv the /tmp file + the
ppg_full <- modify_ppg(ppg_raw)              #   scientificNameAuthorship join load_ppg does
ppg_db   <- format_ppg_for_ts(ppg_full)
ref      <- ts_parse_names(unique(ppg_db$scientificName), tbl_out = TRUE, quiet = TRUE)
mm       <- read_csv("_targets/user/data_raw/pterido_manual_match.csv", show_col_types = FALSE)
```

Run scripts from the project root (`/home/rstudio/ftol/ftol`) so `renv` is picked up — a
`cd scratchpad && Rscript …` in one command fails with "no package called 'targets'".

---

## Failure 1 — `modify_ppg()` has stale workarounds

`modify_ppg()` (in `R/functions.R`) does three kinds of patch, each a candidate for removal
once PPG catches up. Comments like `# will be fixed in ppg v 0.0.0.9008` are the loudest
tells, but check every entry.

- **`filter(taxonID != "wfo-...")`** — drops a duplicate record. Still needed only if the
  name is *still* duplicated in a way that makes matching ambiguous. Check:
  `grep '<name>' /tmp/ppg_new.csv` — if there's now a single accepted record (or the stray
  copy became `unchecked`/`synonym` and `format_ppg_for_ts()` will drop it anyway), delete
  the filter. Confirm the "drop it anyway" case by diffing `ppg_db` with and without the
  filter — if identical, the filter is dead.
- **`dct_add_row(...)`** — adds a taxon missing from PPG. `grep` the new csv for the
  `scientificName`. If it's there now with a compatible `taxonRank` / `taxonomicStatus`,
  delete the `dct_add_row` (leaving it risks a duplicate-name error from dwctaxon).
- **`dct_modify_row(...)`** — forces a status/rank. If the new csv already has that status
  (e.g. `taxonomicStatus == "accepted"`), delete it. If it's *still* `unchecked`, keep it
  and update the stale comment.

Keep entries whose name is still genuinely absent (`grep` returns nothing) or still
`unchecked`. If you add a brand-new `dct_add_row` for a taxon not in PPG, you need the
parent `taxonID` — `awk -F, '$2=="<Genus>" && $4=="genus"' /tmp/ppg_new.csv`.

After editing, run `modify_ppg(load_ppg("0.0.0.9xxx"))` standalone — dwctaxon will error
on the spot if an added row now collides.

### Special case — a name is `nom. ined.` / `nom. inval.` with no valid combination

Sometimes a species turns up (usually via the Failure 6 checkpoint) where NCBI's name
isn't just "missing from PPG" — it's a name PPG can't accept as written because no one has
formally published it in its correct genus yet (NCBI often flags this itself, e.g.
`Terpsichore pacifica (nom. inval.)` / `Lellingeria reunionensis (nom. ined.)`). The tree
usually still makes it obvious which genus it belongs in (BLAST/placement puts it
unambiguously in a genus other than the one in its old name, and it isn't a mis-ID). The
default should be **include it**, not drop it — dropping loses real GenBank data over a
nomenclatural technicality — but *how* to include it depends on whether anyone else has
already put it in that genus, even informally:

- **Someone else has already used a "comb. ined." for it** (check the literature — a
  subsequent paper, checklist, or database entry using the new genus name for this epithet,
  even informally) — adopt that name. PPG itself uses exactly this convention (`grep
  'comb\. ?ined' /tmp/ppg_new.csv`, e.g. `Abacopteris afra (Christ) comb. ined.`,
  `Abrodictyum truncatum (Copel.) comb. ined.`). Follow it exactly: `<New genus> <epithet>
  (<original epithet-publishing author(s)>) comb. ined.` — the parenthesized author is
  whoever published the *epithet* (the basionym author), not whoever used the informal
  combination. In `modify_ppg()`:
  ```r
  dwctaxon::dct_add_row(
    scientificName = "Mycopteris pacifica (Sundue) comb. ined.",
    taxonomicStatus = "accepted",
    taxonRank = "species",
    nomenclaturalStatus = "valid",   # matches PPG's own comb. ined. rows
    parentNameUsageID = "<new genus taxonID>",
    stamp_modified = FALSE
  ) |>
  ```
  Then add manual_match rows in `pterido_manual_match.csv` mapping the NCBI query name(s)
  to that new PPG name — the genus changed, so nothing will fuzzy-match on its own. Add
  **both** forms NCBI carries for these: the full authored name (`Terpsichore pacifica
  Sundue`) and the bracketed species-only form (`Terpsichore pacifica (nom. inval.)`) —
  both can appear as independent queries in `pterido_names_to_inspect`, mirroring how the
  CSV already carries two rows for `Whittieria hengduanensis` (with/without author).
  Verify via `ts_match_names()`/`ts_resolve_names()` as in Failure 2, and confirm the
  standalone `pterido_names_to_inspect` reproduction (Failure 6) comes back 0 rows for
  these queries with `match_type == "manual"`.
- **Nobody has used that combination anywhere, informally or otherwise** — **do not coin
  one ourselves.** "comb. ined." is a real bibliographic claim (someone, somewhere, used
  this name) and FTOL isn't the one to make it. Instead, keep the taxon under its
  *original* name in `modify_ppg()`:
  ```r
  dwctaxon::dct_add_row(
    scientificName = "Terpsichore pacifica Sundue",  # NOT "Mycopteris pacifica … comb. ined."
    taxonomicStatus = "accepted",
    taxonRank = "species",
    stamp_modified = FALSE
  ) |>
  ```
  This needs no `pterido_manual_match.csv` entry — it matches NCBI's own name exactly, so
  it resolves `exact` rather than `manual`. Its accessions will then legitimately break the
  genus it phylogenetically falls into (an old-genus-named tip sitting inside the new
  genus's clade), so add that genus to `make_taxa_exclude_tbl()` (`R/functions.R`, used by
  `non_mono_check`/`verify_non_mono_taxa`) with a one-line comment and issue link — do this
  for the genus the *sequence* falls into (e.g. `Mycopteris`), not the one in the old name.
  File an issue on **github.com/fernphy/ftol** (not PPG — this isn't a PPG problem, it's an
  FTOL workaround) recording the accessions affected, why, and what to undo once a
  combination is published; link it from both the `modify_ppg()` comment and the
  `make_taxa_exclude_tbl()` entry. One issue per species is cleaner than bundling — they'll
  likely be resolved by different, unrelated future publications.
- **Excluding the accessions entirely** (`_targets/user/data_raw/accs_exclude.csv`) is a
  last resort — only when the placement itself is uncertain (not just unnamed), or Joel
  says so for a specific case.

### Related, different case — a valid combination exists but PPG hasn't caught up

Don't confuse the above with a name that *has* been validly, formally published (a real
paper makes the combination) but PPG simply hasn't accepted it into the backbone yet — e.g.
`Zealandia powellii (Baker) Testo & A.R.Field` (Syst. Bot. 44(4): 749, 2019), still a
synonym of `Microsorum powellii` in PPG. No "comb. ined." involved; just add the accepted
name via `dct_add_row`/`dct_modify_row` and sink the old name as its synonym, citing the
publication in the comment. This genuinely is PPG's problem to fix, so file the issue at
**pteridogroup/ppg/issues**, not fernphy/ftol, and link it from the comment.

### Filing a PPG issue: use the taxonomic request form, not a plain issue

PPG issues for a species-level (or below) change go through the **taxonomic request**
template, not a freehand issue — `gh issue create` doesn't render `.yml` issue forms, so
build the body to match what the form would produce:

1. Get the current field list straight from the repo (it can change):
   `gh api repos/pteridogroup/ppg/contents/.github/ISSUE_TEMPLATE/taxonomic-request.yml
   --jq '.content' | base64 -d`. (A change at genus-or-higher instead uses
   `taxonomic-proposal.yml` — same idea, different template; check that one's fields too if
   the request is above species level.)
2. Write the body as `### <field label>` headers in the template's order, one per field,
   each followed by the answer on its own line(s) — that's the markdown GitHub itself
   generates from a filled-in form. Required fields as of this writing: name of taxon
   (with authors — comma-separate if multiple), WFO ID (comma-separate; write "not yet in
   WFO/PPG" for a name that doesn't have one, e.g. a new combination), type of change (pick
   from the dropdown's options — comma-separate if several apply), description of change,
   and the Code of Conduct checkbox (`- [x] I agree to follow the PPG Code of Conduct`).
3. **Always fill "Name of the person(s) making the request" with something like "Written
   by Claude (Anthropic AI) on behalf of `<github username>`"** — never write the request
   as if a human wrote it unprompted.
4. Pass `-l "taxonomic request"` to `gh issue create` (the label the form applies
   automatically isn't added by a plain `gh` call).
5. `gh` may not be installed in this container (`apt-get install -y gh` fixes that), and
   its config may live somewhere other than `$HOME/.config/gh` (e.g. mounted from the host
   at `/home/rstudio/ftol/.config/gh`) — check for an existing `hosts.yml` there before
   assuming a fresh `gh auth login` is needed; point `gh` at it with `GH_CONFIG_DIR=<dir>
   gh auth status` (prefix every `gh` call, since env vars don't persist between Bash
   calls in this harness).

If you (or a past session) already filed a plain-body issue for one of these, close it
(`gh issue close <n> -R pteridogroup/ppg -r "not planned"`, with a comment pointing to the
replacement) and refile properly — don't leave both open.

## Failure 2 — `ts_match_names()`: "manually matched reference names not in reference data"

Every `match` value in `pterido_manual_match.csv` **and** the inline `plastome_manual_match`
list must exist verbatim in `ppg_ref_names$name` (i.e. `unique(ppg_db$scientificName)`,
authorship included). PPG renames things between releases (orthography `-us`→`-os`,
adding `×` hybrid signs, author-string edits) and the old target vanishes. This one error
blocks **both** `plastome_metadata_renamed` and `match_results_raw_round_1`.

Find the offenders:

```r
setdiff(mm$match, ref$name)                                    # CSV
setdiff(<inline plastome_manual_match>$match, ref$name)        # inline tibble
```

For each offender, decide:

- **Query now matches on its own** — test
  `ts_match_names("<query>", unique(ppg_db$scientificName), max_dist = 5,
  match_no_auth = TRUE, match_canon = TRUE, collapse_infra = TRUE, simple = TRUE)`.
  If `match_type` is `exact` / `auto_noauth` / `auto_punct` / `auto_basio*` / `auto_exin-`
  (anything but `auto_fuzzy` / `no_match`), delete the CSV row.
- **Target was renamed** — `grep` the new csv for the taxon, update `match` to the current
  string (e.g. `Lindsaeosoria × flynnii W.H.Wagner` → `× Lindsaeosoria × flynnii W.H.Wagner`).
- **Target became `unchecked`** — see Failure 3.

`ts_match_names()` also refuses duplicate `query` values in the combined list. If you point
a CSV row at a new target, make sure a plastome-inline row for the same query isn't left
pointing elsewhere.

## Failure 3 — a manual-match target is an `unchecked` PPG name

`format_ppg_for_ts()` keeps only `taxonomicStatus %in% c("accepted","synonym")` (plus a few
`nomenclaturalStatus` exclusions), so an `unchecked` name never reaches `ppg_ref_names` and
can't be a `match` target. If you *want* to resolve a query to that name (rather than to
whatever accepted name it should fold into), add its `taxonID` to the `recs_keep` tibble
inside `format_ppg_for_ts()`, with a `# Scientific name` comment like the existing rows.

Get the taxonID: `grep '<name>' /tmp/ppg_new.csv`. Note `format_ppg_for_ts()` builds
`scientificName` as `paste(scientificName, scientificNameAuthorship)` squished — a
no-author PPG record (e.g. `Christella procurrens` with empty authorship) becomes just
`Christella procurrens` in the reference, and the CSV `match` must be written that way.

Verify: after the edit, `"<name incl author>" %in% ts_parse_names(unique(
format_ppg_for_ts(modify_ppg(ppg_raw))$scientificName))$name` is `TRUE`, and
`ts_resolve_names()` returns a non-NA `resolved_name` for it (an `unchecked` name resolves
to itself with `resolved_status == "unchecked"`, which is fine here).

## Failure 4 — inline `plastome_manual_match` is stale / over-stuffed

`resolve_pterido_plastome_names()` has a `plastome_manual_match <- tibble(query=…, match=…)`
that predates most of `pterido_manual_match.csv` and drifts. To find exactly which inline
rows still earn their place: run the plastome match with **no** manual list and list the
`auto_fuzzy` / `no_match` queries; the inline list only needs rows for queries that are
still unresolved *and* not already handled by the CSV.

```r
q  <- unique(<plastome_names_query>$query_name)   # rebuild as in the function, or
                                                 # tar_read a cached copy if fresh enough
bad <- ts_match_names(q, ref, max_dist = 5, match_no_auth = TRUE, match_canon = TRUE,
                      collapse_infra = TRUE, simple = TRUE) |>
  filter(!str_detect(query, " sp\\.$"), str_detect(match_type, "fuzzy|no_match")) |>
  distinct(query) |> pull(query)
setdiff(bad, mm$query)     # → these must stay in the inline list; drop the rest
```

Then confirm the trimmed inline list + CSV gives `resolve_pterido_plastome_names()` zero
`fuzzy|no_match` and zero NA rows (the function's own `assert(not_na, everything())` /
`verify(!any(str_detect(match_type, "fuzzy|no_match")))` enforce this).

## Failure 5 — `create_patel_inclusion_list()`: `assert(not_na, matched_name)` fails

The Patel et al. 2019 Thelypteridaceae supplement uses old spellings; a `case_when` on
`raw_name` in `create_patel_inclusion_list()` rewrites them to current PPG names before
`ts_resolve_names()`. When PPG changes an epithet (seen: `Sphaerostephanos heterocarpus` →
`heterocarpos`, `polycarpus` → `polycarpos`) an entry goes stale or a new one is needed.
Note taxastand's fuzzy match does **not** absorb a trailing `-us`/`-os` change on a short
epithet — it comes back `no_match`, so it must be an explicit `case_when` line.

Diagnose: run `create_patel_inclusion_list(path_to_patel_data = contentid::resolve(
"hash://sha256/233607dc3945dc0f764c44d1171f8bd8bdfe50c4028c9c44e82965e5a5f11fdc",
registries = "local.tsv"), tax_ref = <fresh ppg_db>)`; it errors on the first `no_match`.
Or inspect: `ts_resolve_names(query = <patel_accs$raw_name>, ref_taxonomy = <ppg_db>, …) |>
filter(is.na(matched_name) | str_detect(match_type, "fuzzy"))`. Add/fix the
`raw_name == "<old>" ~ "<current PPG name incl author>"` line; the Patel `raw_name` is the
`accepted_name` column of sheet 3 with underscores replaced by spaces.

## Failure 6 — `ncbi_accepted_names_map`: the human-review checkpoint

This is **not** a bug — `make_ncbi_accepted_names_map(strict = TRUE)` deliberately refuses
to run while `pterido_names_to_inspect` has rows. `inspect_ts_results()` populates it from
every `auto_fuzzy` and `no_match` in `match_results_resolved_all`. Read it:
`tar_read(pterido_names_to_inspect)`. Each `query` is an NCBI name; `grep` it in
`/tmp/ppg_new.csv` and classify:

| What you find in the new ppg.csv | Fix |
|---|---|
| PPG added a `×` hybrid sign; NCBI still has the plain name (`query_match_taxon_agree == TRUE`, `taxonRemarks == "author variant"`) | add `pterido_manual_match.csv` row: `<NCBI name>,<PPG × name>,<note>` |
| Name present but `taxonomicStatus == "unchecked"`, with a clear accepted equivalent | CSV row → the **accepted** name (e.g. `Stegnogramma centrochinensis …` → `Leptogramma centrochinensis Ching ex Y.X.Lin`). Or, to keep the literal name, `recs_keep` (Failure 3). |
| Name genuinely absent from PPG, but the NCBI taxid has real GenBank sequences worth keeping (often a recently-described species; NCBI may tag it `(nom. ined.)`) | `dct_add_row()` in `modify_ppg()` — `taxonomicStatus = "accepted"`, `taxonRank = "species"`, `parentNameUsageID = "<genus taxonID>"`, `stamp_modified = FALSE`. Ask Joel to confirm accepted-vs-synonym if unclear. |
| Absent from PPG, no data worth keeping | add the taxid's accessions to `_targets/user/data_raw/accs_exclude.csv` |

The fuzzy matcher often lists several wrong candidates for one query plus the right `×`
one — resolving the one real target clears all the rows for that query. To confirm you've
cleared it before a full run, reproduce `inspect_ts_results(combined_match_results(
ncbi_names_query = tar_read(ncbi_names_query), r1, r2, r3))` with `r1/r2/r3` built from the
edited CSV + rebuilt `ppg_db`, and check it's 0 rows.

`working_pterido_names_to_inspect.csv` in the repo root is Joel's scratch notes from past
rounds — useful for seeing how similar cases were handled, not read by the pipeline.

---

To actually launch and monitor the fixed-up pipeline (log location, background launch,
stale-lock/zombie recovery, and the `mpcheck_monophy` hang that's since been fixed), switch
to the `run-tar-make` skill.

## Optional — file the underlying PPG problems upstream

Genuine PPG data-quality issues found along the way (duplicate records, `accepted` +
`unchecked` pairs for one name, malformed author strings, missing taxa) are worth filing
so the workaround can be dropped next time — use the taxonomic request form as described
above, not a plain issue. A scan for the common ones:

```r
p <- read_csv("/tmp/ppg_new.csv", show_col_types = FALSE) |>
  mutate(nm = str_squish(paste(scientificName, coalesce(scientificNameAuthorship, ""))))
# same name+author, >1 accepted, or accepted+synonym, or accepted+unchecked:
p |> add_count(nm) |> filter(n > 1) |> group_by(nm) |>
  filter(any(taxonomicStatus == "accepted") &&
         n_distinct(taxonomicStatus) > 1 || sum(taxonomicStatus == "accepted") > 1) |>
  ungroup() |> arrange(nm) |>
  select(taxonID, nm, taxonRank, taxonomicStatus, nomenclaturalStatus, acceptedNameUsageID)
```
