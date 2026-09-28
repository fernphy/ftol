# publish_figshare.R ----
#
# Upload the FigShare-hosted archive files (restez_sql_db.tar.gz,
# taxdmp.zip, README.genbank, README.txt) to the FTOL input data deposit
# (https://doi.org/10.6084/m9.figshare.19474316), replacing the manual SFTP
# download + FigShare web UI upload/delete dance in docs/updating.md
# steps 5-7.
#
# Uses upload_to_figshare_verified() (R/functions.R), which already
# overwrites, so there is no separate "delete old file" step. It wraps
# upload_to_figshare() with checksum verification and retry: a broken pipe
# on a large upload leaves a stale partial file on FigShare that blocks a
# clean retry (cleaned up automatically), and a successful upload can print
# a cosmetic error (github.com/ropenscilabs/deposits/issues/99) that's
# otherwise indistinguishable from a real failure without checking the
# remote checksum directly.
#
# Run this after the targets workflow has finished (same precondition as
# R/snapshot_ftol_data.R -- both just read from the same finished
# _targets store; order between the two doesn't matter).

library(targets)
library(fs)
source("R/packages.R")
source("R/functions.R")

# FigShare deposit ID for FTOL input data
# https://doi.org/10.6084/m9.figshare.19474316
ftol_input_data_deposit_id <- 19474316

# Load targets ----
tar_load(
  c(
    restez_sql_db_archive,
    taxdump_zip_file,
    gb_readme_path,
    input_data_readme
  )
)

# input_data_readme renders under its own filename; FigShare expects
# README.txt (docs/updating.md step 6). tar_render() returns c(output,
# source Rmd) -- take the rendered output, matching the [[1]] pattern in
# R/snapshot_ftol_data.R.
readme_txt <- path(path_temp(), "README.txt")
file_copy(input_data_readme[[1]], readme_txt, overwrite = TRUE)

# Upload files ----
files_to_upload <- c(
  restez_sql_db_archive,
  taxdump_zip_file,
  gb_readme_path,
  readme_txt
)

for (f in files_to_upload) {
  message("Uploading ", f, " to FigShare deposit ", ftol_input_data_deposit_id)
  upload_to_figshare_verified(f, ftol_input_data_deposit_id)
}

file_delete(readme_txt)

message("FigShare publish complete -- all files verified by checksum.")
