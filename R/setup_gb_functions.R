#' Get the index for plants section of GenBank in current GenBank release
#' 
#' Only a single index should be returned; will error if not.
#' 
get_plants_index <- function() {

  temp_file <- tempfile()
  url <- 'https://ftp.ncbi.nlm.nih.gov/genbank/gbrel.txt'
  curl::curl_download(url = url, destfile = temp_file)
  
  downloadable_table <- identify_downloadable_files(temp_file)
  types <- sort(table(downloadable_table[['descripts']]), decreasing = TRUE)
  
  plants_index <- grep("plant|Plant", names(types))
  
  assertthat::assert_that(
    length(plants_index) == 1,
    msg = "Multiple indexes detected for plants"
  )
  
  fs::file_delete(temp_file)
  
  return(plants_index)
}

#' Parse GenBank release downloadable files list
#'
#' Modified from function of same name in restez to accept any path as input
#' 
#' @name identify_downloadable_files
#' @title Identify downloadable files
#' @description Searches through the release notes
#' for a GenBank release to find all listed .seq files.
#' Returns a data.frame for all .seq files and their
#' description.
#' @param flpath Path to release notes downloaded from NCBI FTP server
#' @return data.frame
#' @family private
identify_downloadable_files <- function(flpth) {
  lines <- readLines(con = flpth)
  filesize_section <- filesize <- kill_switch <- descript <-
    descript_section <- FALSE
  filesize_lines <- descript_lines <- NULL
  for (line in lines) {
    if (grepl(pattern = '^[0-9\\.]+\\sFile Descriptions', x = line)) {
      descript_section <- TRUE
      next
    }
    if (grepl(pattern = '^File Size\\s+File Name', x = line)) {
      filesize_section <- TRUE
      next
    }
    if (grepl(pattern = '^[0-9]+\\.\\s', x = line)) {
      descript <- TRUE
    } else {
      descript <- FALSE
    }
    if (grepl(pattern = '^(\\s+)?[0-9]+\\s+gb[a-z]{1,4}[0-9]{1,4}\\.seq{0,1}$',
              x = line)) {
      filesize <- TRUE
    } else {
      filesize <- FALSE
    }
    if (descript_section & descript) {
      descript_lines <- c(descript_lines, line)
    }
    if (filesize_section & filesize) {
      filesize_lines <- c(filesize_lines, line)
      kill_switch <- TRUE
    }
    if (kill_switch & line == '') {
      break
    }
  }
  # break up
  pull <- grepl(pattern = '\\.seq', x = descript_lines)
  seq_files_descripts <- sub('^[0-9]+\\.\\s', '', descript_lines[pull])
  seq_files_descripts <- strsplit(x = seq_files_descripts, split = ' - ')
  seq_files <- unlist(lapply(seq_files_descripts, '[', 1))
  descripts <- unlist(lapply(seq_files_descripts, '[', 2))
  descripts <- sub(pattern = ' sequence entries,', replacement = '',
                   x = descripts)
  descripts <- sub(pattern = ' part [0-9]+\\.', replacement = '',
                   x = descripts)
  filesize_info <- strsplit(x = filesize_lines, split = '\\s')
  filesize_info <- lapply(X = filesize_info, function(x) x[x != ''])
  filesizes <- as.numeric(vapply(X = filesize_info, FUN = '[[', i = 1,
                                 FUN.VALUE = character(1)))
  names(filesizes) <- vapply(X = filesize_info, FUN = '[[', i = 2,
                             FUN.VALUE = character(1))
  # repair truncated names (name of flatfile in some cases got truncated
  # e.g. from "gbpln1000.seq" to "gbpln1000.se")
  truncated_names <- names(filesizes)[grepl("\\.se$", names(filesizes))]
  if (length(truncated_names) > 0) {
    names(filesizes)[grepl("\\.se$", names(filesizes))] <-
      paste0(truncated_names, "q")
  }
  res <- data.frame(seq_files = seq_files, descripts = descripts,
             filesizes = filesizes[seq_files])
  if (any(is.na(res))) {
    warning('Not all file information could be ascertained.')
  }
  res
}

#' Get the list of downloadable plant-division .seq files for the current
#' GenBank release
#'
#' Companion to get_plants_index(): instead of an index into the frequency
#' table of division types, returns the actual filenames for the "Plant"
#' division, for use as the per-file branching input in _targets_gb.R.
#'
#' @param depends Dummy argument to force a {targets} dependency edge; unused
#'
#' @return Character vector of .seq filenames
get_plant_seq_files <- function(depends = NULL) {
  temp_file <- tempfile()
  url <- 'https://ftp.ncbi.nlm.nih.gov/genbank/gbrel.txt'
  curl::curl_download(url = url, destfile = temp_file)

  downloadable_table <- identify_downloadable_files(temp_file)
  fs::file_delete(temp_file)

  plants_type <- grep(
    "plant|Plant", unique(downloadable_table[['descripts']]), value = TRUE
  )

  assertthat::assert_that(
    length(plants_type) == 1,
    msg = "Multiple or no matching descriptions detected for plants"
  )

  downloadable_table[
    downloadable_table[['descripts']] == plants_type, 'seq_files'
  ]
}

#' Optionally cap the list of plant division files, for dry runs
#'
#' @param plant_files Character vector of .seq filenames
#' @param cap Numeric; if not NA, truncate plant_files to this many entries
#'
#' @return Character vector of .seq filenames
cap_plant_files <- function(plant_files, cap = NA) {
  if (is.na(cap)) {
    return(plant_files)
  }
  head(plant_files, cap)
}

#' Assert that a new GenBank release is available
#'
#' @param latest_release Numeric; latest release number on NCBI's FTP server
#' @param current_release Numeric; release number of the currently-installed
#'   local database
#'
#' @return TRUE (invisibly) if a new release is available; errors otherwise
assert_new_gb_release <- function(latest_release, current_release) {
  assertthat::assert_that(
    isTRUE(latest_release > current_release),
    msg = "No new GenBank data available; quitting"
  )
  invisible(TRUE)
}

#' Download one GenBank flatfile, filter to target accessions, delete raw file
#'
#' Streaming replacement for restez::db_download() + restez::db_create():
#' downloads a single division file, uses a fast zgrep pre-check
#' (restez:::search_gz()) to skip full parsing if none of acc_filter could
#' possibly be present, and deletes the raw (compressed) flatfile immediately
#' after parsing regardless of outcome -- so peak disk usage is ~1 file at a
#' time instead of the entire division.
#'
#' @param fl Character; base filename (e.g. "gbpln1.seq"), no .gz extension
#' @param acc_filter Character vector of GenBank accessions to keep
#' @param restez_path Path to use as the restez working directory
#' @param max_tries Maximum download attempts before giving up
#'
#' @return Filtered data.frame of matching records, or NULL if none found
download_and_filter_one_file <- function(
  fl, acc_filter, restez_path, max_tries = 5
) {
  fs::dir_create(restez_path, recurse = TRUE)
  restez::restez_path_set(restez_path)

  tries <- 0
  repeat {
    dl_ok <- tryCatch(
      {
        restez:::file_download(fl, overwrite = FALSE)
        TRUE
      },
      error = function(e) {
        message(sprintf("Download of %s failed: %s", fl, conditionMessage(e)))
        FALSE
      }
    )
    if (isTRUE(dl_ok)) break
    tries <- tries + 1
    if (tries >= max_tries) {
      stop(sprintf("Failed to download %s after %d tries", fl, max_tries))
    }
    Sys.sleep(2^tries)
  }

  gz_path <- file.path(restez:::dwnld_path_get(), paste0(fl, ".gz"))
  on.exit(
    if (file.exists(gz_path)) file.remove(gz_path),
    add = TRUE
  )

  has_match <- restez:::search_gz(acc_filter, gz_path)
  if (!isTRUE(has_match)) {
    return(NULL)
  }

  records <- restez:::flatfile_read(gz_path)
  if (length(records) == 0) {
    return(NULL)
  }

  restez:::gb_df_generate(
    records = records,
    min_length = 0, max_length = NULL,
    acc_filter = acc_filter, invert = FALSE
  )
}

#' Build a fresh GenBank database from streamed, filtered records
#'
#' Combines the (possibly-NULL) per-file results from
#' download_and_filter_one_file() and writes them into a brand new database at
#' restez_path -- any pre-existing database there is deleted first, since this
#' always builds from scratch (no incremental partial-database state to
#' protect; that safety instead comes from {targets} caching each per-file
#' branch independently upstream of this step).
#'
#' @param records_list List of data.frames (or NULLs), one per input file
#' @param restez_path Path to use as the restez working directory
#'
#' @return Path to the resulting database file
build_fresh_gb_db <- function(records_list, restez_path) {
  fs::dir_create(restez_path, recurse = TRUE)
  restez::restez_path_set(restez_path)

  db_path <- restez:::sql_path_get()
  if (file.exists(db_path)) {
    restez::restez_disconnect()
    fs::file_delete(db_path)
  }

  combined <- dplyr::bind_rows(purrr::compact(records_list))

  assertthat::assert_that(
    nrow(combined) > 0,
    msg = "No matching records found across any plant division files"
  )

  restez:::gb_sql_add(df = combined)

  db_path
}

#' Copy a scratch-built GenBank file to its official data_raw location
#'
#' Keeps a single ".bak" copy of the outgoing file before overwriting -- once
#' a GenBank release is superseded, NCBI's FTP server only serves the current
#' release's flatfiles, so a filtered database built from a prior release can
#' never be regenerated later. This is a zero-configuration, always-on safety
#' net; see archive_outgoing_gb_db() for longer-term external archival across
#' more than one release.
#'
#' @param src_path Path to the file in the scratch working directory
#' @param dest_path Final destination path (overwritten if it already exists)
#' @param depends Dummy argument to force a {targets} dependency edge; unused
#'
#' @return dest_path
publish_gb_file <- function(src_path, dest_path, depends = NULL) {
  fs::dir_create(fs::path_dir(dest_path), recurse = TRUE)
  if (fs::file_exists(dest_path)) {
    fs::file_copy(dest_path, paste0(dest_path, ".bak"), overwrite = TRUE)
  }
  fs::file_copy(src_path, dest_path, overwrite = TRUE)
  dest_path
}

#' Archive the outgoing GenBank database to external long-term storage
#'
#' Once a GenBank release is superseded, NCBI's FTP server only serves the
#' current release's flatfiles -- the filtered fern database built from a
#' prior release can never be regenerated later, so this preserves more than
#' the single ".bak" copy publish_gb_file() already keeps locally. Copies the
#' outgoing (about-to-be-replaced) restez files into a release-numbered
#' subdirectory of archive_dir, matching the existing gb_release_<N>/ naming
#' convention already used for manual archives. Safe no-op (with a loud
#' warning, since skipping this does mean permanently losing the ability to
#' regenerate the outgoing release's data beyond the local .bak) if
#' archive_dir is NA, doesn't exist, or the copy fails for any reason --
#' never blocks the pipeline from publishing the new release.
#'
#' @param archive_dir Path to the external archive location, or NA to skip
#' @param current_release Numeric; the OLD release number being replaced (in
#'   restez's internal x10 format, e.g. 2720 for release 272.0)
#' @param data_raw Path to the official data_raw directory
#'
#' @return Path to the archive subdirectory, or NA if skipped
archive_outgoing_gb_db <- function(archive_dir, current_release, data_raw) {
  old_restez_dir <- fs::path(data_raw, "restez")
  if (!fs::dir_exists(old_restez_dir)) {
    # Nothing to archive yet (e.g. very first run)
    return(NA_character_)
  }

  if (is.na(archive_dir) || !fs::dir_exists(archive_dir)) {
    warning(
      "No archive_dir configured/reachable: the outgoing GenBank release's ",
      "filtered database will only be kept as a local .bak copy. NCBI does ",
      "not serve old releases' flatfiles, so this cannot be regenerated ",
      "later if that .bak is ever lost."
    )
    return(NA_character_)
  }

  release_3digit <- round(
    current_release / ifelse(nchar(current_release) == 4, 10, 1)
  )
  dest <- fs::path(archive_dir, sprintf("gb_release_%d", release_3digit))

  tryCatch(
    {
      fs::dir_create(dest, recurse = TRUE)
      fs::dir_copy(old_restez_dir, fs::path(dest, "restez"), overwrite = TRUE)
      old_tar <- fs::path(data_raw, "restez_sql_db.tar.gz")
      if (fs::file_exists(old_tar)) {
        fs::file_copy(old_tar, dest, overwrite = TRUE)
      }
      dest
    },
    error = function(e) {
      warning(sprintf(
        "Failed to archive outgoing GenBank release to %s: %s",
        archive_dir, conditionMessage(e)
      ))
      NA_character_
    }
  )
}

#' Download the GenBank release README into the scratch working directory
#'
#' @param restez_path Path to use as the restez working directory
#'
#' @return Path to the downloaded README.genbank
download_gb_readme <- function(restez_path) {
  fs::dir_create(restez_path, recurse = TRUE)
  dest <- file.path(restez_path, "README.genbank")
  download_with_retry(
    "https://ftp.ncbi.nlm.nih.gov/genbank/README.genbank", dest
  )
  dest
}

#' Write the GenBank release number file into the scratch working directory
#'
#' @param latest_release Numeric release number
#' @param restez_path Path to use as the restez working directory
#'
#' @return Path to the written gb_release.txt
write_gb_release <- function(latest_release, restez_path) {
  fs::dir_create(restez_path, recurse = TRUE)
  restez::restez_path_set(restez_path)
  restez:::gbrelease_log(release = latest_release)
  # restez_path_set() nests its working files one level down, in a "restez"
  # subdirectory of the path given to it -- restez_path_get() reflects that
  file.path(restez::restez_path_get(), "gb_release.txt")
}

#' Archive the published restez database files into a single tar.gz for
#' FigShare
#'
#' @param db_path,release_path,readme_path Published restez file paths
#' @param out_path Destination path for the tar.gz archive
#'
#' @return out_path
archive_restez_db <- function(db_path, release_path, readme_path, out_path) {
  fs::dir_create(fs::path_dir(out_path), recurse = TRUE)
  if (fs::file_exists(out_path)) {
    fs::file_delete(out_path)
  }
  # archive_write_files() silently omits nonexistent input files rather than
  # erroring -- check explicitly so a missing input never produces a
  # silently-incomplete archive
  inputs <- c(db_path, release_path, readme_path)
  assertthat::assert_that(
    all(fs::file_exists(inputs)),
    msg = paste(
      "Cannot build archive, missing input file(s):",
      paste(inputs[!fs::file_exists(inputs)], collapse = ", ")
    )
  )
  archive::archive_write_files(
    archive = out_path,
    files = c(db_path, release_path, readme_path),
    format = "tar",
    filter = "gzip"
  )
  out_path
}

#' Download a fresh NCBI taxdmp.zip, backing up the previous one as .bak
#'
#' @param dest_path Final destination path for taxdmp.zip
#' @param depends Dummy argument to force a {targets} dependency edge; unused
#'
#' @return dest_path
download_taxdmp <- function(dest_path, depends = NULL) {
  fs::dir_create(fs::path_dir(dest_path), recurse = TRUE)
  if (fs::file_exists(dest_path)) {
    fs::file_move(dest_path, paste0(dest_path, ".bak"))
  }
  download_with_retry(
    "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip", dest_path
  )
  dest_path
}

#' Notify that a GenBank download has started
#'
#' @param latest_release Numeric release number, interpolated into the body
#' @param depends Dummy argument to force a {targets} dependency edge; unused
#'
#' @return Invisible NULL
send_gb_start_email <- function(latest_release, depends = NULL) {
  send_gb_email(
    subject = "FTOL download started",
    body_html = glue::glue(
      "FTOL downloading of new GenBank release {latest_release} ",
      "has started on {Sys.time()}"
    )
  )
}

#' Notify that a GenBank download has finished
#'
#' @param latest_release Numeric release number, interpolated into the body
#' @param depends Dummy argument to force a {targets} dependency edge; unused
#'
#' @return Invisible NULL
send_gb_done_email <- function(latest_release, depends = NULL) {
  send_gb_email(
    subject = "FTOL download finished",
    body_html = glue::glue(
      "FTOL downloading of new GenBank release {latest_release} has ",
      "finished on {Sys.time()}. Be sure to upload to FigShare and ",
      "update hash in R/setup.R"
    )
  )
}

#' Send an FTOL GenBank-update notification email
#'
#' Shared auth + send logic, called by send_gb_start_email()/
#' send_gb_done_email(). Takes the subject/body as already-resolved strings
#' (rather than interpolating {targets}-tracked values via glue() directly in
#' the plan) so that dependencies like latest_release are visible to
#' {targets}' static dependency scanner as real function arguments, not hidden
#' inside a string literal it can't see into.
#'
#' @param subject Email subject line
#' @param body_html HTML body content
#'
#' @return Invisible NULL
send_gb_email <- function(subject, body_html) {
  email_draft <-
    gmailr::gm_mime() |>
    gmailr::gm_to("joelnitta@gmail.com") |>
    gmailr::gm_from("pteridogroup.no.reply@gmail.com") |>
    gmailr::gm_subject(subject) |>
    gmailr::gm_html_body(body_html)

  options(gargle_oauth_cache = ".secrets")
  secret_json <- list.files(
    ".secrets",
    pattern = "client_secret.*json", full.names = TRUE
  )
  gmailr::gm_auth_configure(path = secret_json)
  gmailr::gm_oauth_client()
  gmailr::gm_auth("pteridogroup.no.reply@gmail.com")

  gmailr::gm_send_message(email_draft)
  invisible(NULL)
}
