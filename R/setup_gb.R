# Setup local copy of GenBank database

library(restez)
library(magrittr)
source("R/_supp_functions.R")
on.exit(restez_disconnect())

# Download data ----
# Specify location to download GenBank database
restez_path_set("/data_raw")
# Download plant database
db_download(preselection = 1)

# Create database ----

# Specify vector of GenBank accessions: all ferns in GenBank between with seq
# length between 10 to 200000 bases until end date.
# 2022/04/15 is cutoff date for GenBank release 249.
# This date needs to be modified for each new release.
fern_query <- format_all_fern_query(end_date = "2022/04/15")
fern_accs <- gb_fetch_accs(fern_query)
# Also include outgroup accessions
og_accs_tbl <- readr::read_csv("_targets/user/data_raw/plastome_outgroups.csv")
keep_accs <- c(fern_accs, og_accs_tbl$accession)

restez_connect()
db_create(acc_filter = keep_accs, scan = TRUE)

# Copy database to FTOL folder ----
fs::dir_create("_targets/user/data_raw/restez")
fs::dir_copy("/data_raw/restez/sql_db", "_targets/user/data_raw/restez/sql_db")
restez_disconnect()