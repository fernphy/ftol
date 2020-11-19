# Run this script to archive raw data each time the project version is bumped

# Setup ----

# Load functions, packages, and the drake cache
source("_drake.R")

# Prepare metadata ----

# Make tibble of raw data files to copy into the zip archive
raw_data_to_copy <-
  plan %>% 
  filter(format == "file") %>%
  transmute(target, path = as.character(command)) %>%
  mutate(copy = TRUE)

# Make tibble of raw data files *not* to copy into the zip archive
# (large files that won't change; ie, raw sequencing data)

# - load metadata of GoFlag files
loadd(seq_cap_samples, cache = ftol_cache)

raw_data_dont_copy <- 
  tibble(
    path = list.files(
      "data_raw/goflag", 
      recursive = TRUE, 
      pattern = paste(seq_cap_samples, collapse = "|"), 
      full.names = TRUE)
  ) %>%
  mutate(copy = FALSE)

# Make tibble with path of raw data readme to copy into the zip archive
readme <-
  tibble(
    path = list.files("docs", pattern = "raw_data_README.md", full.names = TRUE),
    copy = TRUE
  )

# Combine into single tibble
raw_data_meta <-
  bind_rows(
    raw_data_to_copy,
    raw_data_dont_copy,
    readme
  ) %>%
  mutate(
    file = fs::path_file(path),
    hash = tools::md5sum(path)) %>%
  select(target, file, path, copy, hash)

# Archive the data ----

archive_raw_data(
  version = "ftol_raw_data_v0.0.1", # Change version as needed
  out_path = "working/figshare_raw_data_versions",
  metadata = raw_data_meta)
