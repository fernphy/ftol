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

#' Provide the `contentid` hashes for input data files
#'
#' The hashes in the list should match those of the corresponding data files
#' obtained with `contentid::content_id()` for the dataset at the DOI given.
#'
#' For code to register URLs to hashes, see ./R/register.R
#'
#' @param doi DOI for FTOL input data repository
#' @return List of hashes.
doi_to_hashes <- function(doi) {
  switch(
    doi,
    "10.6084/m9.figshare.19474316.v1" = list(
      accs_exclude = "hash://sha256/e8e215870706a03e74f21a7e117246d77ff1029d91599cda53ec14ea7fbcc1ab", #nolint
      equisetum_subgenera = "hash://sha256/a93ec0663a65d687921af3c412279034786fba769d73408c432bd9b738bd37ad", #nolint
      plastome_outgroups = "hash://sha256/36bf35dbe61c4f133ba5b7112681316b4338480fb6b1299b762f893d5e89c6d1", #nolint
      ppgi_taxonomy_mod = "hash://sha256/d94cf3b3230a4fafaadf76b355a9d989cc1645467aab47934a73cba2920fff3f", #nolint
      target_coding_genes = "hash://sha256/304cd16b67e1d4f180624bed3b683c848c42b6ad5d8250cda1f0425e58831ccf", #nolint
      ref_aln_tar_gz = "hash://sha256/3c37fb9478a8d6d1d5cf12652f01e04c3187db64923be824ca689d924facde18" #nolint
    ),
    "Not a valid DOI"
  )
}

# Load targets ----
tar_load(
  c(
    # paths to input data files
    accs_exclude_path,
    equisteum_subgen_path,
    plastome_outgroups_path,
    ppgi_taxonomy_path,
    target_plastome_genes_path,
    ref_aln_archive,
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
    ftolr_readme
  )
)

# Pre-commit checks ----

# - Check targets status
targets_status <- tar_git_status_targets()

assert_that(
  nrow(targets_status) == 0,
  msg = "One or more outdated targets detected.")

# - Check code status
code_status <- tar_git_status_code()

assert_that(
  nrow(code_status) == 0,
  msg = "Code repo not clean.")

# - Check data
# Expect data hashes to match those registered with contentid for this DOI
# For code to register data, see ./R/register.R
data_doi <- "10.6084/m9.figshare.19474316.v1"
hashes <- doi_to_hashes(data_doi)

assert_that(content_id(accs_exclude_path) == hashes$"accs_exclude")
assert_that(content_id(equisteum_subgen_path) == hashes$equisetum_subgenera)
assert_that(content_id(plastome_outgroups_path) == hashes$plastome_outgroups)
assert_that(content_id(ppgi_taxonomy_path) == hashes$ppgi_taxonomy_mod)
assert_that(
  content_id(target_plastome_genes_path) == hashes$target_coding_genes)
assert_that(content_id(ref_aln_archive) == hashes$ref_aln_tar_gz)

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
  data={data_doi}
  docker={image_tag}")
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
