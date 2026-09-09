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
    send_gb_email(
      subject = "FTOL download started",
      body_html = glue::glue(
        "FTOL downloading of new GenBank release {latest_release} ",
        "has started on {Sys.time()}"
      ),
      depends = plant_files
    )
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

  # Publish to official data_raw location ----
  tar_target(
    restez_db_published,
    publish_gb_file(gb_db_path, path(data_raw, "restez/sql_db")),
    format = "file"
  ),
  tar_target(
    gb_release_published,
    publish_gb_file(
      gb_release_scratch_path, path(data_raw, "restez/gb_release.txt")
    ),
    format = "file"
  ),
  tar_target(
    gb_readme_published,
    publish_gb_file(
      gb_readme_scratch_path, path(data_raw, "restez/README.genbank")
    ),
    format = "file"
  ),

  # Archive for FigShare ----
  tar_target(
    restez_tar_archive,
    archive_restez_db(
      restez_db_published, gb_release_published, gb_readme_published,
      path(data_raw, "restez_sql_db.tar.gz")
    ),
    format = "file"
  ),

  # Download new taxdmp.zip, backing up the previous one ----
  # (done after GenBank data so taxonomic data stays consistent with it)
  tar_target(
    taxdmp_path,
    download_taxdmp(path(data_raw, "taxdmp.zip"), depends = restez_db_published),
    format = "file"
  ),

  # Notify download is finished ----
  done_email = if (send_email_setting) {
    send_gb_email(
      subject = "FTOL download finished",
      body_html = glue::glue(
        "FTOL downloading of new GenBank release {latest_release} has ",
        "finished on {Sys.time()}. Be sure to upload to FigShare and ",
        "update hash in R/setup.R"
      ),
      depends = list(restez_tar_archive, taxdmp_path)
    )
  }
)
