# GenBank ----

#' Format a query string to search for fern sequences from GenBank
#'
#' Returns string to query GenBank for fern (Polypodiopsida) sequences between
#' 10 and 200,000 bp within the given time window
#' 
#' @param start_date String; start date to include sequences
#' @param end_date String; end date to include sequences
#'
#' @return String
format_all_fern_query <- function(start_date = "1980/01/01", end_date) {
  glue::glue(
    'Polypodiopsida[ORGN] AND 10:200000[SLEN] AND ("{start_date}"[PDAT]:"{end_date}"[PDAT]) AND GenBank[filter]' # nolint
  )
}

#' Retrieve GenBank accessions
#'
#' @param query Character vector of length 1; query string to search GenBank
#' @param strict Logical vector of length 1; should an error be issued if
#' the number of accessions retrieved does not match the number of hits
#' from GenBank
#' @return Character vector; accession numbers resulting from query
gb_fetch_accs <- function(query, strict = TRUE) {

  # Conduct search and keep results on server,
  # don't download anything yet
  search_res <- rentrez::entrez_search(
    db = "nuccore",
    term = query,
    use_history = TRUE,
    retmax = 0
  )

  # Make sure something is in there
  assertthat::assert_that(
    search_res$count > 0,
    msg = "Query resulted in no hits"
  )

  # Number of hits NCBI allows us to download at once.
  # This should not need to be changed
  max_hits <- 10000

  # NCBI won't return more than 10,000 results at a time.
  # So download in chunks to account for this
  if (search_res$count > max_hits) {

    # Determine number of chunks
    n_chunks <- search_res$count %/% max_hits

    # Set vector of start values: each chunk
    # will be downloaded starting from that point
    start_vals <- c(0, seq_len(n_chunks) * max_hits)

    # Loop over start values and download up to max_hits for each,
    # then combine
    accessions <- purrr::map_chr(
      start_vals, ~ rentrez::entrez_fetch(
        db = "nuccore",
        web_history = search_res$web_history,
        rettype = "acc",
        retstart = .,
        retmax = max_hits
      )
    ) %>%
      paste(collapse = "")
  } else {
    accessions <- rentrez::entrez_fetch(
      db = "nuccore",
      web_history = search_res$web_history,
      rettype = "acc"
    )
  }

  # NCBI returns accessions as single string, so split into vector
  accessions <- strsplit(x = accessions, split = "\\n")[[1]]
  accessions <- sub(pattern = "\\.[0-9]+", replacement = "", x = accessions)

if (strict) {
  n_accs <- length(accessions)
  assertthat::assert_that(
    search_res$count == n_accs,
    msg = glue::glue(
      "Number of accessions ({n_accs}) not equal to number of GenBank hits ({search_res$count})" # nolint
    )
  )
}

  accessions
}

#' Get a tibble of sequence lengths and descriptions from a local GenBank
#' database
#'
#' @param accessions Character vector; GenBank accession numbers.
#' @param restez_path Path to the GenBank database prepared with {restez}.
#' @param na_rm Logical vector of length 1; should any NA values in accessions
#' be removed?
#'
#' @return Tibble with three columns, `accession`, `def`, and `slen`
#'
gb_slen_get <- function(accessions, restez_path = "/data_raw", na_rm = TRUE) {

  accessions <- unique(accessions)
  assertthat::assert_that(assertthat::is.flag(na_rm))
  if (na_rm) {
    accessions <- accessions[!is.na(accessions)]
  }

  # Connect to restez database
  suppressMessages(restez::restez_path_set(restez_path))
  suppressMessages(restez::restez_connect())
  # Disconnect from restez database when done
  on.exit(restez::restez_disconnect())

  # Extract sequence definitions from database
  defs <- restez::gb_definition_get(accessions)

  # Extract sequence lengths from database
  seq_lens <-
    restez::gb_sequence_get(accessions) %>% nchar()

  # Combine into single tibble
  def_tib <-
    tibble(
      accession = names(defs),
      def = defs
    ) %>%
    assert(is_uniq, accession)

  seq_len_tib <-
    tibble(
      accession = names(seq_lens),
      slen = seq_lens
    ) %>%
    assert(is_uniq, accession)

  def_tib %>%
    left_join(seq_len_tib, by = "accession") %>%
    assert(is_uniq, accession) %>%
    assert(not_na, everything())
}

# Data inspection ----

#' Get useful information about accessions for a species in FTOL
#'
#' @param species_select Species name.
#' @param sanger_accessions_selection Dataframe of accessions used in FTOL.
#' @param ncbi_accepted_names_map Dataframe mapping NCBI taxid to resolved
#'   name from pteridocat.
#' @param match_results_resolved_all Dataframe with results of resolving
#'   NCBI taxonomic names to pteridocat.
#'
#' @return Tibble
#'
get_acc_info <- function(
  species_select,
  sanger_accessions_selection,
  ncbi_accepted_names_map,
  match_results_resolved_all) {
  sanger_accessions_selection %>%
    filter(species == species_select) %>%
    mutate(across(contains("seq_len"), as.character)) %>%
    left_join(
      unique(select(ncbi_accepted_names_map, species, taxid, resolved_name)),
      by = "species"
    ) %>%
    left_join(
      select(match_results_resolved_all, -taxid), by = "resolved_name"
    ) %>%
    pivot_longer(everything()) %>%
    mutate(value = as.character(value)) %>%
    filter(!is.na(value)) %>%
    filter(value != "0")
  }

#' Lookup name data in NCBI taxonomic database for a
#' species in FTOL
#'
#' @param species_select Species name.
#' @param sanger_accessions_selection Dataframe of accessions used in FTOL.
#' @param ncbi_accepted_names_map Dataframe mapping NCBI taxid to resolved
#'   name from pteridocat.
#' @param ncbi_names_full NCBI taxonomic names database for ferns
#'
#' @return Tibble
get_ncbi_name <- function(
  species_select,
  sanger_accessions_selection,
  ncbi_accepted_names_map,
  ncbi_names_full) {
  sanger_accessions_selection %>%
    filter(species == species_select) %>%
    select(-contains("seq_len")) %>%
    left_join(
      select(ncbi_accepted_names_map, species, taxid, resolved_name)
    ) %>%
    select(taxid) %>%
    inner_join(ncbi_names_full)
}
