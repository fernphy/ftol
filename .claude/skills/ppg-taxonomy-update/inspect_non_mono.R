# inspect_non_mono.R -- triage a non_mono_check failure (fast-tree monophyly gate).
#
# verify_non_mono_taxa() stops the pipeline when the fast IQ-TREE tree has a
# non-monophyletic taxon that is neither in non_mono_notes.csv nor in
# make_taxa_exclude_tbl(). Most are fast-tree artifacts; a few are real
# mis-IDs the all-by-all BLAST missed; some are naming conflicts (e.g. a PPG
# version change). This prints, per unhandled taxon, what is needed to tell
# them apart. Reads the _targets store, writes nothing.
#
# For each taxon: tips / MRCA clade size / intruder tips (grouped), how many of
# the taxon's tips are NEW vs the previous release's tree (ftol_data/, so renames
# from a PPG bump show up as new + removed tips), and the "drop test": remove the
# listed outlier tips and see whether the taxon (and the ones sharing its
# intruders) become monophyletic. A flag that vanishes when 1-2 tips are dropped
# is usually one misplaced/renamed sample, not a bad taxon.
#
# Usage (repo root; Docker on the bare host):
#   docker run --rm -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) -v $PWD:/wd:ro \
#     -w /wd joelnitta/ftol:latest Rscript .claude/skills/ppg-taxonomy-update/inspect_non_mono.R
#   ... inspect_non_mono.R all     # also show taxa already handled
suppressPackageStartupMessages({library(targets); library(tidyverse); library(ape)})
source("R/functions.R")
show_all <- "all" %in% commandArgs(TRUE)
fm   <- tar_read(fast_monophy_by_clade); nn <- tar_read(non_mono_notes); ex <- make_taxa_exclude_tbl()
tree <- tar_read(sanger_fast_tree_rooted); samp <- tar_read(sanger_sampling)
prev <- read.tree("ftol_data/ftol_sanger_ml.tre")$tip.label   # previous release
new  <- setdiff(tree$tip.label, prev); gone <- setdiff(prev, tree$tip.label)
cat("tips now:", length(tree$tip.label), "| previous release:", length(prev), "| new:", length(new), "| removed:", length(gone), "\n")
cat("NEW tips:", paste(new, collapse = "; "), "\nREMOVED tips:", paste(gone, collapse = "; "), "\n\n")
bad <- filter(fm, monophyly == "No")
if (!show_all) bad <- bad |> anti_join(nn, by = "taxon") |> anti_join(ex, by = "taxon")
cat("taxa to inspect:", nrow(bad), "\n")
levels_ <- c("order", "suborder", "family", "subfamily", "genus", "subgenus")
clade_of <- function(tr, tips) { m <- getMRCA(tr, tips); tr$tip.label[phangorn::Descendants(tr, m, "tips")[[1]]] }
for (tx in bad$taxon) {
  l <- levels_[map_lgl(levels_, ~ tx %in% samp[[.x]])][1]; if (is.na(l)) next
  tips <- intersect(samp$species[samp[[l]] == tx & !is.na(samp[[l]])], tree$tip.label)
  clade <- clade_of(tree, tips); intr <- setdiff(clade, tips)
  out <- str_split(fm$outliers[fm$taxon == tx], ",\\s*")[[1]]; out <- out[out %in% tips]
  grp <- samp |> filter(species %in% intr) |> count(across(any_of(setdiff(levels_, c("order", "suborder"))[1:3])), name = "n") |> arrange(desc(n))
  cat(sprintf("\n===== %s (%s): %d tips, clade %d, intruders %d | new tips in taxon: %d | outliers: %s\n", tx, l, length(tips), length(clade), length(intr), sum(tips %in% new), if (length(out)) paste(out, collapse = ", ") else "-"))
  if (nrow(grp)) { cat("  intruders by group (top):\n"); print(as.data.frame(head(grp, 6)), right = FALSE, row.names = FALSE) }
  if (length(out)) { t2 <- drop.tip(tree, out); tp2 <- setdiff(tips, out); c2 <- clade_of(t2, tp2)
    cat(sprintf("  drop test (remove %d outlier tip(s)): clade %d -> %d, intruders %d -> %d  %s\n", length(out), length(clade), length(c2), length(intr), length(setdiff(c2, tp2)), if (length(setdiff(c2, tp2)) == 0) "=> MONOPHYLETIC without them" else "=> still not")) }
}
