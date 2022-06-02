# GenBank ----

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
gb_slen_get <- function(
  accessions, restez_path = "_targets/user/data_raw", na_rm = TRUE) {

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
