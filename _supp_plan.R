# _supp_plan.R
#
# Workflow to run supplemental analyses other than building the Fern Tree of
# Life (FTOL)

library(targets)
library(tarchetypes)
source("R/packages.R")
source("R/functions.R")
source("R/_supp_functions.R")

Sys.setenv(TAR_PROJECT = "supp_plan")

# Load targets from main workflow
# - Cutoff for sequences to include
date_cutoff <- tar_read(date_cutoff, store = "_targets")
# - Vector of GenBank accessions in local copy of GenBank database
accs_in_local_db <- tar_read(accs_in_local_db, store = "_targets")

tar_plan(
  # Compare seq lengths ----
  # Goal is to plot difference in distribution of seq lengths
  # between whole plastome and Sanger seqs
  #
  # Format query to search for all fern sequences on GenBank
  all_ferns_query = format_all_fern_query(end_date = date_cutoff),
  # Obtain accession numbers
  all_ferns_accs = gb_fetch_accs(all_ferns_query),
  # Limit scope to fern accessions in local GenBank db
  fern_accs = all_ferns_accs[all_ferns_accs %in% accs_in_local_db$accession],
  # Assemble tibble of sequence descriptions and lengths
  fern_slen_tib = gb_slen_get(fern_accs)
)