# Setup local copy of GenBank database for ferns
library(fs)
library(restez)
library(gmailr)
library(assertthat)

source("R/setup_gb_functions.R")

# Should email notifications be sent? ----
send_email <- TRUE

# Compare latest and current GenBank release ----
# - get latest release number on FTP server
latest_release <- restez:::latest_genbank_release() |>
  as.numeric()

# - get current downloaded release number
# (assumes we've already done this once)
current_gb_release_file <- "_targets/user/data_raw/restez/gb_release.txt"

current_release <- readLines(current_gb_release_file) |>
  as.numeric()

# - run checks
invisible(
  assert_that(
    file_exists(current_gb_release_file),
    msg = "Cannot find current GenBank release file; quitting"
  )
)

invisible(
  assert_that(
    isTRUE(latest_release > current_release),
    msg = "No new GenBank data available; quitting"
  )
)

invisible(
  assert_that(
    !file_exists("scratch/dl_running.txt"),
    msg = "Download currently in progress; quitting"
  )
)

# Get plants index
plants_index <- get_plants_index()

# Prepare temporary download folder ----
# DELETES OLD DATA (flatfiles)
if (dir_exists("scratch")) {
  dir_delete("scratch")
}
dir_create("scratch")
# presence of scratch/dl_running.txt means that download is in progress
writeLines(as.character(latest_release), "scratch/dl_running.txt")

# Send email to indicate that download has started
# (setup credentials with R/setup_email.R)
if (send_email) {
  # Draft email
  email_draft <-
    gm_mime() |>
    gm_to("joelnitta@gmail.com") |>
    gm_from("pteridogroup.no.reply@gmail.com") |>
    gm_subject("FTOL download started") |>
    gm_html_body(
      glue::glue(
        "FTOL downloading of new GenBank release {latest_release} \\
        has started on {Sys.time()}"
      )
    )

  # Authenticate email server
  options(gargle_oauth_cache = ".secrets")
  secret_json <- list.files(
    ".secrets",
    pattern = "client_secret.*json", full.names = TRUE
  )
  gm_auth_configure(path = secret_json)
  gm_oauth_client()
  gm_auth("pteridogroup.no.reply@gmail.com")

  # Send email
  gm_send_message(email_draft)
}

# Download data ----
# Specify location to download GenBank flatfiles and create database
restez_path_set("scratch")
# Download plant flatfiles
# this includes >900 files totaling >600 gb, so may get interrupted during dl
# use max_tries to automatically restart
db_download(preselection = plants_index, overwrite = FALSE, max_tries = 1000)

# Create database ----
# Specify vector of GenBank accessions:
# all ferns in GenBank between with seq ength between 10 to 200000 bases
# (longest fern plastome < 200000 bp)
fern_accs <- ncbi_acc_get("Polypodiopsida[ORGN] AND 10:200000[SLEN]")
# Also include outgroup accessions
og_accs_tbl <- readr::read_csv("_targets/user/data_raw/plastome_outgroups.csv")
keep_accs <- unique(c(fern_accs, og_accs_tbl$accession))
# Create db
db_create(acc_filter = keep_accs, scan = TRUE)
# Download genbank README
download.file(
  "https://ftp.ncbi.nlm.nih.gov/genbank/README.genbank",
  "scratch/restez/README.genbank"
)

# Copy database (not flatfiles) to FTOL folder ----
# will overwrite old data
if (dir_exists("_targets/user/data_raw/restez")) {
  dir_delete("_targets/user/data_raw/restez")
}
dir_create("_targets/user/data_raw/restez")
file_copy(
  "scratch/restez/sql_db",
  "_targets/user/data_raw/restez/sql_db",
  overwrite = TRUE
)
file_copy(
  "scratch/restez/gb_release.txt",
  "_targets/user/data_raw/restez/gb_release.txt",
  overwrite = TRUE
)
file_copy(
  "scratch/restez/README.genbank",
  "_targets/user/data_raw/restez/README.genbank",
  overwrite = TRUE
)

# Also compress to tar archive for figshare
restez_tar <- "_targets/user/data_raw/restez_sql_db.tar.gz"
if (file_exists(restez_tar)) {
  file_delete(restez_tar)
}
archive::archive_write_files(
  archive = restez_tar,
  files = c(
    "scratch/restez/sql_db",
    "scratch/restez/gb_release.txt",
    "scratch/restez/README.genbank"
  ),
  format = "tar",
  filter = "gzip"
)

# Cleanup ----

# Download done, so delete "running" file
# (scratch dir including all GenBank flatfiles remains, handle manually)
if (file_exists("scratch/dl_running.txt")) {
  file_delete("scratch/dl_running.txt")
}

if (send_email) {
  # Draft email
  email_draft <-
    gm_mime() |>
    gm_to("joelnitta@gmail.com") |>
    gm_from("pteridogroup.no.reply@gmail.com") |>
    gm_subject("FTOL download started") |>
    gm_html_body(
      glue::glue(
        "FTOL downloading of new GenBank release {latest_release} has \\
        finished on {Sys.time()}. Be sure to upload to FigShare and update \\
        hash in R/setup.R"
      )
    )

  # Authenticate email server
  options(gargle_oauth_cache = ".secrets")
  secret_json <- list.files(
    ".secrets",
    pattern = "client_secret.*json", full.names = TRUE
  )
  gm_auth_configure(path = secret_json)
  gm_oauth_client()
  gm_auth("pteridogroup.no.reply@gmail.com")

  # Send email
  gm_send_message(email_draft)
}
