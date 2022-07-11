# Setup local copy of GenBank database

library(restez)
library(magrittr)
source("R/functions.R")
on.exit(restez_disconnect())

# Download data ----
# Specify location to download GenBank database
restez_path_set("scratch")

# Download plant database
# Connection may get dropped, resulting in an error.
# Repeat this in a while() loop until it completes
tries <- 0
while (TRUE) {
  x <- try(db_download(preselection = 1))
  if (inherits(x, "try-error")) {
    cat("ERROR: ", x, "\n")
    tries <- tries + 1
    message(paste("Trying again, attempt number", tries))
    Sys.sleep(10)
   } else {
    break
   }
}

# Create database ----

# Specify vector of GenBank accessions:
# all ferns in GenBank between with seq ength between 10 to 200000 bases
# (longest fern plastome < 200000 bp)
fern_accs <- gb_fetch_accs("Polypodiopsida[ORGN] AND 10:200000[SLEN]")
# Also include outgroup accessions
og_accs_tbl <- readr::read_csv("_targets/user/data_raw/plastome_outgroups.csv")
keep_accs <- unique(c(fern_accs, og_accs_tbl$accession))

restez_connect()

db_create(acc_filter = keep_accs, scan = TRUE)

# Copy database to FTOL folder ----
# will overwrite old data
fs::dir_delete("_targets/user/data_raw/restez")
fs::dir_create("_targets/user/data_raw/restez")
fs::file_copy(
  "scratch/restez/sql_db",
  "_targets/user/data_raw/restez/sql_db",
  overwrite = TRUE)
fs::file_copy(
  "scratch/restez/gb_release.txt",
  "_targets/user/data_raw/restez/gb_release.txt",
  overwrite = TRUE)

# Also compress to tar archive for figshare
archive::archive_write_files(
  archive = "_targets/user/data_raw/restez_sql_db.tar.gz",
  files = c("scratch/restez/sql_db", "scratch/restez/gb_release.txt"),
  format = "tar",
  filter = "gzip"
)

restez_disconnect()