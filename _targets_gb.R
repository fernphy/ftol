# Load packages and functions
library(targets)
library(tarchetypes)
library(crew)
library(restez)
library(assertthat)
# NB: gmailr is NOT attached on purpose. gmailr 2.0 exports a defunct
# `message()` that would mask base::message() (which R/setup_gb_functions.R
# and R/functions.R call directly), turning any logged download failure into
# a spurious "message() is defunct" error. All gmailr use is via `gmailr::`.
source("R/packages.R")
source("R/functions.R")
source("R/setup_gb_functions.R")

Sys.setenv(TAR_PROJECT = "gb_download")

# Specify location of raw data (final, official location read by _targets.R)
# Overridable for dry runs so test output never touches real project data
data_raw <- Sys.getenv("GB_DL_DATA_RAW", "_targets/user/data_raw")

# Working directory for streaming download + filter (kept small: each plant
# division file is deleted immediately after it's scanned/parsed)
# Overridable for dry runs, same reason as above
scratch_dir <- Sys.getenv("GB_DL_SCRATCH", "scratch")

# Should email notifications be sent? (disable for dry runs)
send_email_setting <- as.logical(Sys.getenv("GB_DL_SEND_EMAIL", "TRUE"))

# Optional cap on number of plant division files to process, for dry runs
# validating the pipeline end-to-end without the full multi-day download
plant_files_cap <- as.numeric(Sys.getenv("GB_DL_FILE_CAP", NA))

# External location for long-term archival: after each release is published,
# its filtered database is copied to <dir>/gb_release_<N>/ (one snapshot per
# release). NCBI only serves the current release's flatfiles, so a superseded
# release's filtered database can never be regenerated -- this is the only
# copy of past releases. Unset/NA skips archival (with a loud warning).
archive_dir_setting <- Sys.getenv("GB_DL_ARCHIVE_DIR", NA)

# Set options:
# - Modest local parallelization: empirically validated safe against NCBI's
#   file server (6 concurrent downloads, no throttling/errors), kept well
#   below the main pipeline's 20 workers to be a good FTP citizen
# - Continue past errored branches so one flaky file doesn't block the rest;
#   simply re-invoke tar_make() later to retry any that failed
tar_option_set(
  error = "continue",
  controller = crew_controller_local(workers = 4)
)

tar_plan(
  # Check for new release ----
  # Always re-check on every invocation, since this watches external state
  tar_target(
    latest_release,
    gb_release_number(),
    cue = tar_cue(mode = "always")
  ),
  tar_file_read(
    current_release,
    path(data_raw, "restez/gb_release.txt"),
    as.numeric(readLines(!!.x))
  ),
  # Errors (halting the whole pipeline) if no new release is available
  release_check = assert_new_gb_release(latest_release, current_release),

  # Get list of plant-division files to process ----
  # (optionally capped to a handful of files for dry runs, see plant_files_cap)
  plant_files = cap_plant_files(
    get_plant_seq_files(depends = release_check), plant_files_cap
  ),

  # Determine which accessions to keep ----
  tar_file_read(
    outgroup_accs,
    path(data_raw, "plastome_outgroups.csv"),
    read_csv(!!.x)
  ),
  fern_accs = ncbi_acc_get("Polypodiopsida[ORGN] AND 10:200000[SLEN]"),
  keep_accs = unique(c(fern_accs, outgroup_accs$accession)),

  # Notify download has started (fires once per new release; cached on resume)
  start_email = if (send_email_setting) {
    send_gb_start_email(latest_release, depends = plant_files)
  },

  # Stream download + filter, one plant-division file at a time ----
  tar_target(
    plant_file_records,
    download_and_filter_one_file(plant_files, keep_accs, scratch_dir),
    pattern = map(plant_files)
  ),

  # Build fresh database from filtered records ----
  # deployment = "main": this target depends on the whole plant_file_records
  # pattern (thousands of branches). Dispatching that to a crew worker with
  # the default retrieval = "main" makes the main process spin at 100% CPU
  # leaking memory instead of transferring the data -- the target never
  # completes. Running it in the main process sidesteps the transfer; the
  # build itself is a ~5-second, <1 GB operation.
  tar_target(
    gb_db_path,
    build_fresh_gb_db(plant_file_records, scratch_dir),
    format = "file",
    deployment = "main"
  ),

  # Download GenBank README and write release number to scratch ----
  # These stay in scratch until the very end. The official gb_release.txt is
  # the gate release_check reads, so it must be the LAST official file to
  # change -- a crash any time before that leaves release_check seeing the
  # old release, and a resumed run correctly picks the work back up.
  tar_target(
    gb_readme_scratch_path,
    download_gb_readme(scratch_dir),
    format = "file"
  ),
  tar_target(
    gb_release_scratch_path,
    write_gb_release(latest_release, scratch_dir),
    format = "file"
  ),

  # Bundle for FigShare ----
  # Built from the scratch copies (not the published ones) so it needn't wait
  # on publishing -- it only needs the new release's data to exist somewhere
  tar_target(
    restez_tar_archive,
    archive_restez_db(
      gb_db_path, gb_release_scratch_path, gb_readme_scratch_path,
      path(data_raw, "restez_sql_db.tar.gz")
    ),
    format = "file"
  ),

  # Download new taxdmp.zip ----
  # (after the GenBank data, so taxonomic data stays consistent with it)
  tar_target(
    taxdmp_path,
    download_taxdmp(path(data_raw, "taxdmp.zip"), depends = gb_db_path),
    format = "file"
  ),

  # Publish database + README to the official data_raw location ----
  # gb_release.txt is deliberately NOT published here -- it goes last (below).
  tar_target(
    restez_db_published,
    publish_gb_file(
      gb_db_path, path(data_raw, "restez/sql_db"),
      depends = list(restez_tar_archive, taxdmp_path)
    ),
    format = "file"
  ),
  tar_target(
    gb_readme_published,
    publish_gb_file(
      gb_readme_scratch_path, path(data_raw, "restez/README.genbank"),
      depends = list(restez_tar_archive, taxdmp_path)
    ),
    format = "file"
  ),

  # Archive this release to external long-term storage ----
  # After the db/README/bundle exist but before the gb_release.txt gate flips,
  # so a failed archive still leaves release_check seeing the old release.
  # Never errors: returns NA (with a warning, surfaced in the done email) if
  # the archive dir is unset/unreachable.
  archive_result = archive_gb_db(
    archive_dir_setting, latest_release,
    restez_db_published, gb_readme_published, restez_tar_archive
  ),

  # Flip the release gate LAST ----
  tar_target(
    gb_release_published,
    publish_gb_file(
      gb_release_scratch_path, path(data_raw, "restez/gb_release.txt"),
      depends = list(
        restez_db_published, gb_readme_published, archive_result,
        restez_tar_archive, taxdmp_path
      )
    ),
    format = "file"
  ),

  # Notify download is finished ----
  done_email = if (send_email_setting) {
    send_gb_done_email(
      latest_release, archive_result, depends = gb_release_published
    )
  }
)
