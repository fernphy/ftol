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
  filesizes <- as.integer(vapply(X = filesize_info, FUN = '[[', i = 1,
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
