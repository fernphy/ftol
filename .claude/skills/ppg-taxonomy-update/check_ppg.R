# check_ppg.R -- pre-flight for a PPG version, BEFORE launching tar_make().
#
# Rebuilds ppg_full -> ppg_db -> reference from the *current* R/functions.R for
# the given PPG tag and runs every name-resolution gate we know about, using
# the upstream targets already in the _targets store (tar_read): stale
# modify_ppg() workarounds, manual-match targets missing from the reference,
# the plastome gate, the sanger-side pterido_names_to_inspect checkpoint, the
# higher-taxonomy build (ppg_tl / ppgi_taxonomy) and create_patel_inclusion_list()
# (the skill's Failure 5). Writes nothing. Each stage prints [ OK ] or [FAIL]
# and a failing stage does not abort the later ones.
#
# Usage (from the repo root; needs the pipeline's R packages, i.e. Docker on
# the bare host, or the dev container's shell):
#   docker run --rm -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) \
#     -v $PWD:/wd:ro -w /wd joelnitta/ftol:latest \
#     Rscript .claude/skills/ppg-taxonomy-update/check_ppg.R 2.0.0
# (The plastome/sanger stages need ncbi_names_query and the plastome metadata
# targets, so run it after a pipeline run has got that far, e.g. after a
# failure at plastome_metadata_renamed.)

ver <- commandArgs(TRUE)[1]
suppressPackageStartupMessages({library(targets); library(tidyverse); library(taxastand); library(assertr); library(readr); library(dwctaxon); library(taxlist)})
source("R/functions.R")
verify <- function(data, ...) data   # plastome gate disabled: we want to SEE offenders
step <- function(label, expr) {
  tryCatch({ r <- expr; cat(sprintf("  [ OK ] %s\n", label)); r },
    error = function(e) { cat(sprintf("  [FAIL] %s: %s\n", label, substr(gsub("\n", " ", conditionMessage(e)), 1, 300))); NULL })
}
cat("######## PPG", ver, "########\n")
ppg_raw  <- step("load_ppg", load_ppg(ver))
ppg_full <- step("modify_ppg (stale workarounds collide?)", modify_ppg(ppg_raw))
if (is.null(ppg_full)) quit(save = "no")
ppg_db   <- step("format_ppg_for_ts", format_ppg_for_ts(ppg_full))
ref      <- ts_parse_names(unique(ppg_db$scientificName), tbl_out = TRUE, quiet = TRUE)
cat("  ref names:", nrow(ref), " ppg_db rows:", nrow(ppg_db), "\n")
mm <- read_csv("_targets/user/data_raw/pterido_manual_match.csv", show_col_types = FALSE)
miss <- setdiff(mm$match, ref$name)
cat("  manual-match TARGETS missing from reference:", length(miss), "\n"); if (length(miss)) print(head(miss, 15))
# --- plastome
args <- list(plastome_ncbi_names_raw = tar_read(plastome_ncbi_names_raw), plastome_metadata_raw = tar_read(plastome_metadata_raw),
             plastome_outgroups = tar_read(plastome_outgroups), ref_names_parsed = ref, ref_names_data = ppg_db, manual_matches = mm)
src <- readLines("R/functions.R"); st <- grep("^resolve_pterido_plastome_names <- ", src)
en <- st + which(src[(st + 1):length(src)] == "}")[1]; fn <- src[st:en]; cut <- grep("### Combine results ###", fn)
f <- eval(parse(text = paste(c(fn[seq_len(cut - 1)], "  match_results_plastome_resolved", "}"), collapse = "\n")))
pl <- step("plastome resolution", do.call(f, args))
if (!is.null(pl)) { bad <- pl |> filter(str_detect(match_type, "fuzzy|no_match")); cat("  plastome names:", nrow(pl), " fuzzy/no_match:", nrow(bad), "\n"); if (nrow(bad)) print(as.data.frame(select(bad, query, match_type, matched_name)), right = FALSE) }
# --- sanger checkpoint
q <- tar_read(ncbi_names_query); vk <- define_varieties_to_keep()
m <- function(query) ts_match_names(query = query, reference = ref, manual_match = mm, max_dist = 5, match_no_auth = TRUE, match_canon = TRUE, collapse_infra = TRUE, collapse_infra_exclude = vk)
sg <- step("sanger rounds + inspect", {
  q1 <- select_ncbi_names_round_1(q); r1 <- ts_resolve_names(m(q1$scientific_name), ppg_db)
  q2 <- select_ncbi_names_round_2(r1, q); r2 <- ts_resolve_names(m(q2$scientific_name), ppg_db)
  q3 <- select_ncbi_names_round_3(r1, r2, q); r3 <- ts_resolve_names(m(q3$species), ppg_db)
  inspect_ts_results(combined_match_results(ncbi_names_query = q, r1, r2, r3))
})
if (!is.null(sg)) { cat("  pterido_names_to_inspect rows:", nrow(sg), "\n"); if (nrow(sg)) print(as.data.frame(select(sg, query)), right = FALSE) }
# --- higher taxonomy (the ppg_tl / ppgi_taxonomy targets)
tl <- step("dwc_to_tl(ppg_full)", dwc_to_tl(ppg_full))
if (!is.null(tl)) {
  tdf <- step("taxlist_to_df()", taxlist_to_df(tl))
  if (!is.null(tdf)) cat("  ppgi_taxonomy rows:", nrow(tdf), "\n")
}
# --- Patel et al. 2019 inclusion list (Failure 5: stale case_when name repairs)
pat <- step("create_patel_inclusion_list()", create_patel_inclusion_list(
  path_to_patel_data = contentid::resolve(
    "hash://sha256/233607dc3945dc0f764c44d1171f8bd8bdfe50c4028c9c44e82965e5a5f11fdc",
    registries = "local.tsv"),
  tax_ref = ppg_db))
if (!is.null(pat)) cat("  patel inclusion rows:", nrow(pat), "\n")
