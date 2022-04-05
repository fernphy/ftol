# snapshot_ftol_data.R ----
#
# Update ftol_data (**output** of ftol workflow), commit,
# and push to github
#
# This should be run after the targets workflow is done and all
# code changes have been committed.

library(gert)
library(gittargets)
library(assertthat)
library(targets)
library(fs)

# Define helper function ----

#' Write out the CC0 license to a temporary file
#' @param path Path to the file with the CC0 license
#' @return Path where the CC0 license was written
write_cc0 <- function(path) {
  # Load CC0 license from fernphy/ferncal repo
  cc0 <- readr::read_lines("https://raw.githubusercontent.com/fernphy/ferncal/main/LICENSE") # nolint
  # Check it is (probably) what we expect
  assertthat::assert_that(
    cc0[[3]] == "CC0 1.0 Universal",
    msg = "Probably not CC0 license"
  )
  readr::write_lines(cc0, path)
  path
}

# Pre-commit checks ----

# Check targets status
targets_status <- tar_git_status_targets()

assert_that(
  nrow(targets_status) == 0,
  msg = "One or more outdated targets detected.")

# Check code status
code_status <- tar_git_status_code()

assert_that(
  nrow(code_status) == 0,
  msg = "Code repo not clean.")

# Pull from origin ----
git_pull(repo = "ftol_data")

# Add files ----

# Load files to commit
tar_load(c(
  # accessions
  acc_table_long_ftolr, acc_table_wide_ftolr,
  # taxonomy
  sanger_sampling_ftolr,
  # trees
  plastome_tree_ftolr,
  sanger_ml_tree_ftolr, sanger_ml_tree_dated_ftolr,
  sanger_con_tree_ftolr, sanger_con_tree_dated_ftolr,
  # alignments
  sanger_alignment_ftolr, plastome_alignment_ftolr,
  plastome_parts_table_ftolr, sanger_parts_table_ftolr,
  # fossils
  con_fossil_calibration_tips_ftolr, ml_fossil_calibration_tips_ftolr,
  # README
  ftolr_readme
))

# Write out CC0 license
write_cc0("ftol_data/LICENSE")

# Add files to commit
ftol_data_files <-
  path_file(
    c(
    # accessions
    acc_table_long_ftolr, acc_table_wide_ftolr,
    # taxonomy
    sanger_sampling_ftolr,
    # trees
    plastome_tree_ftolr,
    sanger_ml_tree_ftolr, sanger_ml_tree_dated_ftolr,
    sanger_con_tree_ftolr, sanger_con_tree_dated_ftolr,
    # alignments
    sanger_alignment_ftolr, plastome_alignment_ftolr,
    plastome_parts_table_ftolr, sanger_parts_table_ftolr,
    # fossils
    con_fossil_calibration_tips_ftolr, ml_fossil_calibration_tips_ftolr,
    # README
    ftolr_readme[[1]],
    "LICENSE"
    )
  )

added <- git_add(
  files = ftol_data_files,
  repo = "ftol_data"
)

# Commit changes and push ----

if (nrow(added) > 0) {
  # Format message:
  # commit hash of code repo
  # plus comment of code repo
  msg <- glue::glue("code={git_commit_info()$id}
  comment={git_commit_info()$message}")
  # Commit
  git_commit(
    repo = "ftol_data",
    message = msg
  )
  # Push
  git_push(remote = "origin", repo = "ftol_data")
} else {
  print("No changes to add; nothing committed or pushed")
}
