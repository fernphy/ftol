# Setup local copy of GenBank database for ferns
library(fs)
library(restez)
library(blastula)

# Compare latest and current GenBank release ----
# - get latest release number on FTP server
latest_release <- restez:::latest_genbank_release() |>
  as.numeric()

# - get current downloaded release number
# (assumes we've already done this once)
current_gb_release_file <- "_targets/user/data_raw/restez/gb_release.txt"

assertthat::assert_that(
  file_exists(current_gb_release_file),
  msg = "Cannot find current GenBank release file; quitting"
)

current_release <- readLines(current_gb_release_file) |>
  as.numeric()

# quit if latest release is not newer than current release
if (!isTRUE(latest_release > current_release)) {
  message("No new GenBank data available; quitting")
  quit(save = "no")
}

# quit if a download is currently running
if (fs::file_exists("scratch/dl_running.txt")) {
  message("Download currently in progress; quitting")
  quit(save = "no")
}

# Prepare temporary download folder ----
# DELETES OLD DATA (flatfiles)
if (dir_exists("scratch")) {
  dir_delete("scratch")
}
dir_create("scratch")
# presence of scratch/dl_running.txt means that download is in progress
writeLines(latest_release, "scratch/dl_running.txt")

# Send email to indicate that download has started
# (setup credentials with R/setup_email.R)
email <- compose_email(
  glue::glue(
    "FTOL downloading of new GenBank release {latest_release} \\
    has started on {Sys.time()}")
)

smtp_send(
    email,
    from = "jnitta.no.reply@gmail.com",
    to = "joelnitta@gmail.com",
    subject = "FTOL download started",
    credentials = creds_file(file = "gmail_creds")
  )

# Download data ----
# Specify location to download GenBank flatfiles and create database
restez_path_set("scratch")
# Download plant flatfiles
# this includes >900 files totaling >600 gb, so may get interrupted during dl
# use max_tries to automatically restart
db_download(preselection = 1, overwrite = FALSE, max_tries = 1000)

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
  overwrite = TRUE)
file_copy(
  "scratch/restez/gb_release.txt",
  "_targets/user/data_raw/restez/gb_release.txt",
  overwrite = TRUE)
file_copy(
  "scratch/restez/README.genbank",
  "_targets/user/data_raw/restez/README.genbank",
  overwrite = TRUE)

# Also compress to tar archive for figshare
archive::archive_write_files(
  archive = "_targets/user/data_raw/restez_sql_db.tar.gz",
  files = c(
    "scratch/restez/sql_db",
    "scratch/restez/gb_release.txt",
    "scratch/restez/README.genbank"),
  format = "tar",
  filter = "gzip"
)

# Cleanup ----

# Download done, so delete "running" file
if (fs::file_exists("scratch/dl_running.txt")) {
  fs::file_delete("scratch/dl_running.txt")
}

# Send email to indicate that download has finished
email <- compose_email(
  glue::glue(
    "FTOL downloading of new GenBank release {latest_release} \\
    has finished on {Sys.time()}")
)

smtp_send(
    email,
    from = "jnitta.no.reply@gmail.com",
    to = "joelnitta@gmail.com",
    subject = "FTOL download finished",
    credentials = creds_file(file = "gmail_creds")
  )