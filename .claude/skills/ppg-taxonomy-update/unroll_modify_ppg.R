# unroll_modify_ppg.R -- run modify_ppg() one patch at a time on a PPG tag.
#
# modify_ppg() is one long pipe of filter()/dct_add_row()/dct_modify_row()
# patches; when one fails the error ("tax_dat must be of class data.frame",
# "taxonID in new data must be different from that in existing data", ...)
# doesn't say which. This unrolls the pipe, skips (and reports) failing
# patches and flags NO-OPs, i.e. workarounds PPG has since fixed that can be
# deleted. Writes nothing.
#
# Usage (from the repo root, same Docker/dev-container note as check_ppg.R):
#   Rscript .claude/skills/ppg-taxonomy-update/unroll_modify_ppg.R 2.0.0

suppressPackageStartupMessages({library(tidyverse); library(dwctaxon)})
source("R/functions.R")
ver <- commandArgs(TRUE)[1]
raw <- load_ppg(ver)
b <- body(modify_ppg); e <- b[[length(b)]]
steps <- list()
while (is.call(e) && length(e) >= 2 && is.call(e[[2]]) && !identical(e[[1]], as.name("{"))) { steps <- c(list(e), steps); e <- e[[2]] }
first <- e; steps <- c(list(first), steps)   # innermost: ppg_raw |> filter(...)
cat("steps:", length(steps), "\n")
res <- raw; failed <- list()
for (i in seq_along(steps)) {
  st <- steps[[i]]
  # replace the piped input (first arg) with `res`
  st[[2]] <- quote(res)
  if (i == 1) {  # innermost has ppg_raw symbol as first arg already replaced
  }
  out <- tryCatch(eval(st), error = function(err) structure(list(msg = conditionMessage(err)), class = "step_err"))
  label <- substr(paste(deparse(st[-2]), collapse = " "), 1, 150)
  if (inherits(out, "step_err")) { cat(sprintf("[%2d] ERROR: %s\n     %s\n", i, out$msg, label)); failed[[length(failed) + 1]] <- i }
  else { same <- isTRUE(all.equal(out, res, check.attributes = FALSE)); cat(sprintf("[%2d] %s  %s\n", i, if (same) "NO-OP (dead on this version)" else "changes data", label)); res <- out }
}
cat("failed steps:", paste(unlist(failed), collapse = ","), "\n")
