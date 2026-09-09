# Load packages and functions
library(targets)
library(tarchetypes)
library(crew)
library(restez)
library(gmailr)
library(assertthat)
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

# Optional external location for long-term archival of the outgoing release's
# database before it's overwritten (NCBI only serves the current release's
# flatfiles, so once overwritten, an old release's filtered database can
# never be regenerated). A single local .bak copy is always kept regardless
# of this setting; set this to also copy it somewhere with more headroom for
# multiple past releases (e.g. an external drive), unset/NA to skip.
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
    as.numeric(restez:::latest_genbank_release()),
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
  tar_target(
    gb_db_path,
    build_fresh_gb_db(plant_file_records, scratch_dir),
    format = "file"
  ),

  # Download GenBank README and write release number ----
  # (kept in scratch only for now -- NOT published to the official data_raw
  # location yet. The official gb_release.txt is what release_check reads to
  # decide whether a new release is available, so it must stay untouched
  # until every other step below has actually finished; otherwise a crash
  # partway through would make a resumed run wrongly conclude there's nothing
  # left to do, while the official database was never actually updated.)
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

  # Archive the outgoing release before any official files are overwritten ----
  archive_result = archive_outgoing_gb_db(
    archive_dir_setting, current_release, data_raw
  ),

  # Archive for FigShare ----
  # Built from the scratch copies (not the official published ones) so this
  # doesn't need to wait on publishing -- it only needs the new release's
  # data to exist somewhere, not to already be "official"
  tar_target(
    restez_tar_archive,
    archive_restez_db(
      gb_db_path, gb_release_scratch_path, gb_readme_scratch_path,
      path(data_raw, "restez_sql_db.tar.gz")
    ),
    format = "file"
  ),

  # Download new taxdmp.zip, backing up the previous one ----
  # (done after GenBank data so taxonomic data stays consistent with it)
  tar_target(
    taxdmp_path,
    download_taxdmp(path(data_raw, "taxdmp.zip"), depends = gb_db_path),
    format = "file"
  ),

  # Publish to official data_raw location ----
  # The true final step: only runs once the archive and taxdmp are confirmed
  # done, so gb_release.txt (the release_check gate) is the very last thing
  # to change -- a crash at any point before this leaves the official
  # location, and therefore release_check, untouched
  tar_target(
    restez_db_published,
    publish_gb_file(
      gb_db_path, path(data_raw, "restez/sql_db"),
      depends = list(archive_result, restez_tar_archive, taxdmp_path)
    ),
    format = "file"
  ),
  tar_target(
    gb_readme_published,
    publish_gb_file(
      gb_readme_scratch_path, path(data_raw, "restez/README.genbank"),
      depends = list(archive_result, restez_tar_archive, taxdmp_path)
    ),
    format = "file"
  ),
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
    send_gb_done_email(latest_release, depends = gb_release_published)
  }
)
