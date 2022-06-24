# snapshot_ftol_data.R ----
#
# Update ftol_data (**output** of ftol workflow), commit,
# and push to github remote https://github.com/fernphy/ftol_data
#
# This should be run after the targets workflow is done and all
# code changes have been committed.

library(gert)
library(gittargets)
library(assertthat)
library(targets)
library(fs)
library(contentid)

# Define helper functions ----

#' Write out the CC0 license
#' @param path Path to write CC0 license
#' @return Path where the CC0 license was written
write_cc0 <- function(path) {
  cc0_file <- contentid::resolve("hash://sha256/6d489af6292662d9e36d34ce49423784984a5f6e41d7b58f49b01264df59fa03") # nolint
  fs::file_move(cc0_file, path)
  path
}

# Load targets ----
tar_load(
  c(
    # paths to input data files
    ref_aln_archive,
    restez_sql_db_archive,
    # docker tag
    image_tag,
    # output data to commit
    # - accessions
    acc_table_long_ftolr, acc_table_wide_ftolr,
    # - taxonomy
    sanger_sampling_ftolr,
    # - trees
    plastome_tree_ftolr,
    sanger_ml_tree_ftolr, sanger_ml_tree_dated_ftolr,
    sanger_con_tree_ftolr, sanger_con_tree_dated_ftolr,
    # - alignments
    sanger_alignment_ftolr, plastome_alignment_ftolr,
    plastome_parts_table_ftolr, sanger_parts_table_ftolr,
    # - fossils
    con_fossil_calibration_tips_ftolr, ml_fossil_calibration_tips_ftolr,
    # - README
    ftol_data_readme
  )
)

# Pre-commit checks ----

# - Check targets status (no targets should be outdated)
targets_status <- tar_git_status_targets()

assert_that(
  nrow(targets_status) == 0,
  msg = "One or more outdated targets detected.")

# - Check code status (no uncommited changes should be present)
code_status <- tar_git_status_code()

# - Check data
# Only need to check data files that are archived outside of this repo
# (ie, on figshare https://doi.org/10.6084/m9.figshare.19474316)
# If this fails, run content_id(ref_aln_archive) to obtain new hash.
assert_that(
  content_id(ref_aln_archive) == "hash://sha256/388b53201a8626d4b41851e716505e7904d24ee3730de25310cb82cd3a1e6e71", # nolint
  msg = glue::glue("Contents of {ref_aln_archive} have changed; update hash")
  )
assert_that(
  content_id(restez_sql_db_archive) == "hash://sha256/8059a845c6570eeffb6fe08c29e178a9dc223ab6f929a1b6c6b374e160f21410", # nolint
  msg = glue::glue(
    "Contents of {restez_sql_db_archive} have changed; update hash")
  )

# - Check renv status (lock file is synchronized with library)
assert_that(
  isTRUE(renv::status()$synchronized),
  msg = "Renv lock file not in sync with library"
)

# Pull from origin ----
git_pull(repo = "ftol_data")

# Add files ----

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
    ftol_data_readme[[1]],
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
  data={data_doi}
  docker={image_tag}")
  # Commit
  git_commit(
    repo = "ftol_data",
    message = msg
  )
  # Push
  git_push(remote = "origin", repo = "ftol_data", verbose = TRUE)
} else {
  print("No changes to add; nothing committed or pushed")
}
