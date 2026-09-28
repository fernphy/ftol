# publish_figshare.R ----
#
# Upload the FigShare-hosted archive files (restez_sql_db.tar.gz,
# taxdmp.zip, README.genbank, README.txt) to the FTOL input data deposit
# (https://doi.org/10.6084/m9.figshare.19474316), replacing the manual SFTP
# download + FigShare web UI upload/delete dance in docs/updating.md
# steps 5-7.
#
# Uses upload_to_figshare() (R/functions.R), which already overwrites, so
# there is no separate "delete old file" step.
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
# README.txt (docs/updating.md step 6)
readme_txt <- path(path_temp(), "README.txt")
file_copy(input_data_readme, readme_txt, overwrite = TRUE)

# Upload files ----
files_to_upload <- c(
  restez_sql_db_archive,
  taxdump_zip_file,
  gb_readme_path,
  readme_txt
)

for (f in files_to_upload) {
  message("Uploading ", f, " to FigShare deposit ", ftol_input_data_deposit_id)
  upload_to_figshare(f, ftol_input_data_deposit_id)
}

file_delete(readme_txt)

message("FigShare publish complete.")
