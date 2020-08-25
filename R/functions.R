# Download GenBank seqs ----

#' Helper function to parse character vector from GenBank record into dataframe
#'
#' @param text Character vector from reading in a GenBank record
#'
#' @return Dataframe
#' 
#' @examples
#' gb_text <- c("ACCESSION   MK697585",                                             
#' "REFERENCE   1  (bases 1 to 817)",                                  
#' "  TITLE     Hybridization rates in Dryopteris carthusiana complex",
#' "REFERENCE   2  (bases 1 to 817)",                                  
#' "  TITLE     Direct Submission")  
#' tidy_gb_text(gb_text)
tidy_gb_text <- function(text) {
  
  if(length(text) == 0) return(tibble())
  if(is.null(text)) return(tibble())
  
  text %>%
    stringr::str_trim("left") %>%
    # Split text on the first instance of >1 space
    stringr::str_split("  +", n = 2) %>%
    # Convert to tibble in wide format
    purrr::map_df(~tibble::tibble(var = .[[1]], value = .[[2]])) %>%
    dplyr::mutate(var = janitor::make_clean_names(var)) %>%
    tidyr::pivot_wider(values_from = "value", names_from = "var")
}

#' Fetch reference data (title of reference where sequence was published)
#' for a GenBank query
#'
#' @param query GenBank query
#'
#' @return Tibble
#' @export
#'
#' @examples
#' fetch_genbank_refs("KY241392")
#' fetch_genbank_refs("Crepidomanes minutum AND rbcL[Gene]")
fetch_genbank_refs <- function(query) {
  
  # Get list of GenBank IDs (GIs)
  uid <- reutils::esearch(term = query, db = "nucleotide", usehistory = TRUE)
  
  # Exract number of hits and print
  num_hits <- reutils::content(uid, as = "text") %>% str_match("<eSearchResult><Count>([:digit:]+)<\\/Count>") %>% magrittr::extract(,2)
  print(glue("Found {num_hits} sequences (UIDs)"))
  
  # Download complete GenBank record for each and write it to a temporary file
  temp_dir <- tempdir()
  temp_file <- fs::path(temp_dir, "gb_records.txt")
  
  reutils::efetch(uid, "nucleotide", rettype = "gb", retmode = "text", outfile = temp_file)
  
  
  # Read in full GenBank records text
  gbrecs <- readr::read_lines(temp_file) %>% 
    # Replace problematic slashes (that indicate a break between records) with "RECORD_BREAK"
    stringr::str_replace_all("^//$", "RECORD_BREAK") %>% 
    # Keep only relevant lines in each record
    magrittr::extract(str_detect(., "ACCESSION|TITLE|REFERENCE|RECORD_BREAK"))
  
  # Make a vector of record groups for splitting the records into a list
  record_groups <- cumsum(gbrecs == "RECORD_BREAK")
  
  # Split the records into a list of records
  split(gbrecs, record_groups) %>%
    # Don't need the "RECORD_BREAK" string anymore
    purrr::map(~magrittr::extract(., stringr::str_detect(., "RECORD_BREAK", negate = TRUE))) %>%
    # Tidy the list
    purrr::map_df(tidy_gb_text)
  
}

#' Extract a translated amino acid sequence from a single entry
#' in a genbank flatfile
#'
#' @param gb_entry String (character vector of length 1);
#' single entry from a genbank flatfile
#' @param gene Name of gene
#'
#' @return amino acid sequence as a character vector, named
#' for the species + voucher
#' 
extract_sequence <- function (gb_entry, gene) {
  
  # Extract start and end of target gene
  gene_range_list <-
    gb_entry %>%
    paste(sep = "") %>%
    str_remove_all("\n") %>%
    str_remove_all('\"') %>%
    str_match("FEATURES(.*)ORIGIN") %>%
    magrittr::extract(,1) %>%
    # Match strings like 'gene complement(<1..10) /gene=rbcL' 
    # (we want the gene and the range it contains)
    # use negative look-ahead to match middle part NOT containing the word "/gene"
    # use '?' for non-greedy match, that will stop on the first instance of "/gene"
    str_match_all("gene +(?!/gene).*?/gene=[:alnum:]+") %>%
    unlist
  
  # Make sure target gene is detected
  assertthat::assert_that(
    any(str_detect(gene_range_list, regex(gene, ignore_case = TRUE))),
    msg = "Gene not detected"
  )
  
  # Subset to only the target gene
  gene_range <-
    gene_range_list %>%
    magrittr::extract(str_detect(., regex(gene, ignore_case = TRUE))) %>% 
    str_split(" +") %>%
    unlist %>%
    magrittr::extract(!str_detect(., "gene")) %>%
    magrittr::extract(str_detect(., "\\d")) %>%
    str_match_all("\\d+") %>%
    unlist() %>%
    parse_number %>%
    sort()
  
  # Make sure that worked correctly
  assertthat::assert_that(
    length(unique(gene_range)) <= 2,
    msg = "Duplicate copies of gene detected")
  assertthat::assert_that(
    length(gene_range) > 1,
    msg = "Full range of gene not detected")
  assertthat::assert_that(is.numeric(gene_range))
  assertthat::assert_that(!anyNA(gene_range))
  assertthat::assert_that(gene_range[1] <= gene_range[2])
  assertthat::assert_that(gene_range[2] >= gene_range[1])
  
  # Extract sequence, subset to target gene
  sequence <- 
    gb_entry %>%
    paste(sep = "") %>%
    str_remove_all("\n") %>%
    str_remove_all('\"') %>%
    str_match('ORIGIN(.+)$') %>%
    magrittr::extract(,2) %>%
    str_remove_all(" ") %>%
    str_remove_all("[0-9]") %>%
    substr(gene_range[1], gene_range[2])
  
  # Extract accession
  accession <-
    gb_entry %>%
    paste(sep = "") %>%
    str_match('ACCESSION(.+)\n') %>%
    magrittr::extract(,2) %>%
    # In very rare cases, may have multiple values for accession,
    # separated by space. If so, take the first one.
    str_trim(side = "both") %>%
    str_split(" ") %>%
    purrr::pluck(1,1)
  
  set_names(sequence, accession)
}

#' Parse a genbank file and extract all the DNA sequences for a given gene
#'
#' @param gbff_path Path to genbank flat file
#' @param gene Name of gene
#'
#' @return List of class "AAbin"
#' 
parse_dna_from_flatfile <- function (gbff_path, gene) {
  
  # Read-in flat file
  readr::read_file(gbff_path) %>%
    # '\\' is delimiter between entries
    stringr::str_split("\n\\/\\/\n") %>%
    unlist %>%
    # Drop the last item, as it is just an empty line (after the last '\\')
    magrittr::extract(-length(.)) %>%
    # Extract DNA sequences from each entry
    purrr::map2(gene, extract_sequence) %>%
    # Name them as the accession
    purrr::set_names(map_chr(., names)) %>%
    # Convert to ape format
    # - each sequence needs to be a character vector with each letter as an element
    map(stringr::str_split, "") %>%
    map(ape::as.DNAbin) %>%
    jntools::flatten_DNA_list()
  
}

#' Download a set of fern sequences for a given gene
#' 
#'
#' @param gene Name of gene
#' @param start_date Earliest date to download
#' @param end_date Most recent date to download
#'
#' @return List of class DNAbin
#' fetch_fern_gene("rbcL", start_date = "2018/01/01", end_date = "2018/01/10")
fetch_fern_gene <- function(gene, start_date = "1980/01/01", end_date) {
  
  assertthat::assert_that(assertthat::is.string(gene))
  
  assertthat::assert_that(assertthat::is.string(end_date))
  
  # Format GenBank query: all ferns matching the gene name.
  # Assume that we only want single genes or small sets of genes, not entire plastome.
  # Set upper limit to 7000 bp (we will fetch plastomes >7000 bp separately).
  query <- glue('{gene}[Gene] AND Polypodiopsida[ORGN] AND 1:7000[SLEN] AND ("{start_date}"[PDAT]:"{end_date}"[PDAT])')
  
  # Get list of GenBank IDs (GIs)
  uid <- reutils::esearch(term = query, db = "nucleotide", usehistory = TRUE)
  
  # Extract number of hits and print
  num_hits <- reutils::content(uid, as = "text") %>% str_match("<eSearchResult><Count>([:digit:]+)<\\/Count>") %>% magrittr::extract(,2)
  print(glue("Found {num_hits} sequences (UIDs)"))
  
  # Download complete GenBank record for each and write it to a temporary file
  temp_dir <- tempdir()
  temp_file <- fs::path(temp_dir, "gb_records.txt")
  
  reutils::efetch(uid, "nucleotide", rettype = "gb", retmode = "text", outfile = temp_file)
  
  # Parse flatfile
  parse_dna_from_flatfile(temp_file, gene)
  
}

#' Download a set of fern sequence metadata for a given gene
#'
#' @param gene Name of gene
#' @param start_date Earliest date to download
#' @param end_date Most recent date to download
#'
#' @return Datarame
#' fetch_fern_metadata("rbcL", start_date = "2018/01/01", end_date = "2018/02/01")
fetch_fern_metadata <- function(gene, start_date = "1980/01/01", end_date) {
  
  assertthat::assert_that(assertthat::is.string(gene))
  
  assertthat::assert_that(assertthat::is.string(end_date))
  
  # Format GenBank query: all ferns matching the gene name and date range.
  # Assume that we only want single genes or small sets of genes, not entire plastome.
  # Set upper limit to 7000 bp (we will fetch plastomes >7000 bp separately).
  query <- glue('{gene}[Gene] AND Polypodiopsida[ORGN] AND 1:7000[SLEN] AND ("{start_date}"[PDAT]:"{end_date}"[PDAT])')
  
  # Fetch standard metadata
  metadata <- gbfetch::fetch_metadata(query) %>%
    # Parse out data contained in the "subtype" and "subname" columns
    tidy_genbank_metadata() %>%
    # Remove weird brackets from species names
    mutate(species = str_remove_all(species, "\\[") %>% str_remove_all("\\]"))
  
  # Also fetch reference data (publication to cite for the sequence).
  # If the title is only "Direct Submission" consider this to be NA
  # (no real pub associated with the sequence, can't match on this)
  ref_data <- fetch_genbank_refs(query) %>%
    transmute(
      accession, 
      publication = str_replace_all(title, "Direct Submission", NA_character_))
  
  # Combine the two
  left_join(metadata, ref_data, by = "accession")
  
}

#' Combine GenBank DNA sequences and associated metadata into single dataframe
#' 
#' Also adds column "length" for the length of the DNA sequence
#'
#' @param seqs DNA sequences (list of class DNAbin)
#' @param metadata Data associated with the sequences. Must include
#' columns "accession" and "species".
#' "accession" must match names of DNA sequences exactly.
#' @param ... Additional arguments not used by this function but meant for 
#' tracking with drake
#'
#' @return Tibble
#'
#' @examples
#' fasta <- fetch_seqs("Crepidomanes minutum psba")
#' fasta_data <- fetch_metadata("Crepidomanes minutum psba")
#' join_genbank_fasta_with_meta(fasta, fasta_data)
join_genbank_fasta_with_meta <- function (seqs, metadata, ...) {
  
  # Check input
  assertthat::assert_that(is.list(seqs))
  assertthat::assert_that(inherits(seqs, "DNAbin"))
  assertthat::assert_that(is.data.frame(metadata))
  
  # Convert sequences to tibble and calculate length
  seqs_tibble <- tibble(seq = as.character(seqs), accession = names(seqs)) %>%
    mutate(length = map_dbl(seq, length))
  
  rm(seqs)
  
  # Check that all sequences are in metata. There should be zero missing, so
  # error if not.
  missing_from_metadata <- anti_join(seqs_tibble, metadata, by = "accession")
  assertthat::assert_that(
    nrow(missing_from_metadata) == 0,
    msg = "Not all sequences in metadata")
  
  # Check for all metadata with sequences. Some seqs may have gotten dropped when
  # extracting gene regions from GenBank flatfiles. Warn if not so.
  missing_from_seqs <- anti_join(metadata, seqs_tibble, by = "accession")
  message <- assertthat::validate_that(
    nrow(missing_from_seqs) == 0,
    msg = warning(glue(
      "The following accessions are in the metadata 
     but missing from sequences, and will be dropped: 
     {paste(missing_from_seqs$accession, collapse = ', ')}")))
  
  # Combine sequences and metadata into single tibble
  # this will make sorting and selecting easier.
  dplyr::inner_join(metadata, seqs_tibble, by = "accession")
  
}

#' Parse species names from GenBank and resolve to World Ferns list automatically
#'
#' @param combined_metadata Tibble including sequences and associated
#' metadata. Sequences in column 'seq' as list-vector.
#' @param col_plants Tibble with taxonomic reference to use
#' for resolving names.
#'
#' @return List including:
#'   - names_resolved_auto: Names that were successfully resolved automatically
#'   to Catalog of Life (i.e., World Ferns)
#'   - names_resolved_to_other_sources: Names that were successfully resolved 
#'   automatically a source other than Catalog of Life
#'   - names_with_mult_syns: Names with multiple synonyms that require
#'   manual inspection and selection
#'   - total_failures: Names that could not be resolved to any of the above
#'   - non_valid_queries: Names that were not valid for querying (not to identified
#'   to species, hybrid, etc.)
#'
#' @examples
#' loadd(col_plants, cache = plastid_cache)
#' fasta <- gbfetch::fetch_sequences("Crepidomanes minutum rbcl")
#' fasta_data <- gbfetch::fetch_metadata("Crepidomanes minutum rbcl")
#' fasta_combined <- combine_rbcl_fasta_with_meta(fasta, fasta_data)
#' resolve_genbank_names_auto(fasta_combined, col_plants)
resolve_genbank_names_auto <- function (combined_metadata, col_plants) {
  
  ### Attempt to resolve names to World Ferns
  # (note that GenBank names only include species without author)
  name_resolution_results <- taxastand::resolve_fern_names(
    names = unique(combined_metadata$species),
    col_plants = col_plants,
    exact_match_by = "species",
    resolve_to = "scientific_name")
  
  # Exclude non-valid names (hybrids, names not identified to species or of
  # uncertain identity) from results
  name_resolution_results_valid_query <-
    name_resolution_results %>%
    filter_at(vars(contains("exclude")), all_vars(. == FALSE)) %>%
    select(-contains("exclude")) %>%
    filter(str_detect(query, " sp. ", negate = TRUE)) %>%
    filter(str_detect(query, "aff.", negate = TRUE)) %>%
    filter(str_detect(query, "cf.", negate = TRUE))
  
  # *add to tally at end
  non_valid_queries <- anti_join(
    name_resolution_results, 
    name_resolution_results_valid_query,
    by = "query")
  
  ### Separate failures and successfully resolved names
  failures <-
    name_resolution_results_valid_query %>%
    filter(is.na(scientificName))
  
  # *add to tally at end
  names_resolved_auto <-
    name_resolution_results_valid_query %>%
    filter(!is.na(scientificName))
  
  ### Make list of species with multiple synonyms. 
  # This will be written out and used to choose names manually.
  
  # - first get CoL database ID
  col_source_id <-
    taxize::gnr_datasources() %>% 
    filter(
      str_detect(title, "Catalogue of Life")
    ) %>%
    pull(id)
  
  # - run query
  names_with_mult_syns <-
    failures %>% 
    filter(fail_reason == "Matches multiple distinct sci names") %>%
    pull(query) %>%
    # re-run gnr_resolve against CoL to get synonyms
    taxize::gnr_resolve(preferred_data_sources = col_source_id, fields = "all") %>% 
    transmute(
      query = user_supplied_name, 
      data_source_title, # To track database version
      imported_at, # To track database version
      matched_name, 
      synonyms = current_name_string,
      name_resolved_manual = NA_character_,
      name_resolved_manual_source = NA_character_
    )
  
  ### Make list of names without any matches in CoL that
  ### can be resolved to GBIF, IPNI, or Tropicos.
  
  # - get source IDs for querying names that lacked any match in CoL.
  # Use GBIF, IPNI, or Tropicos (preferred in that order).
  sources_select <- taxize::gnr_datasources() %>% 
    filter(
      str_detect(title, "GBIF|The International Plant Names Index|Tropicos")
    ) %>%
    pull(id)
  
  # - run query on anything that didn't get resolved to a name with multiple synonyms
  names_resolved_to_other_sources_raw <-
    failures %>%
    anti_join(names_with_mult_syns, by = "query") %>% 
    pull(query) %>%
    taxize::gnr_resolve(preferred_data_sources = sources_select, fields = "all") %>%
    group_by(user_supplied_name) %>%
    add_count(data_source_title)
  
  # - categorize the results into those that had
  # multiple hits within a data source vs. those that didn't
  if(any(names_resolved_to_other_sources_raw$n > 1)) {
    names_resolved_to_other_sources_mult_hits <-
      names_resolved_to_other_sources_raw %>%
      filter(n > 1) %>%
      transmute(
        query = user_supplied_name,
        data_source_title, # To track database version
        imported_at, # To track database version
        matched_name, 
        synonyms = current_name_string,
        name_resolved_manual = NA_character_,
        name_resolved_manual_source = NA_character_
      ) } else {
        names_resolved_to_other_sources_mult_hits <- tibble()
      }
  
  # Combine results of names with mult. synonyms
  # *add to tally at end
  names_with_mult_syns <-
    bind_rows(names_with_mult_syns, names_resolved_to_other_sources_mult_hits)
  
  # Separate out names resolved to a single name in a non-CoL data source
  # *add to tally at end
  names_resolved_to_other_sources <-
    names_resolved_to_other_sources_raw %>%
    anti_join(names_with_mult_syns, by = c(user_supplied_name = "query")) %>%
    filter(n == 1) %>%
    mutate(
      data_source_title = factor(
        data_source_title, 
        levels = c(
          "GBIF Backbone Taxonomy", 
          "The International Plant Names Index",
          "Tropicos - Missouri Botanical Garden")
      )) %>%
    arrange(desc(score), data_source_title) %>%
    slice(1) %>%
    ungroup %>%
    transmute(
      query = user_supplied_name,
      data_source_title, # To track database version
      imported_at, # To track database version
      matched_name, 
      synonyms = ifelse("current_name_string" %in% colnames(.), current_name_string, NA_character_),
      name_resolved_manual = NA_character_,
      name_resolved_manual_source = NA_character_
    )
  
  # Give up on anything left.
  # *add to tally at end
  total_failures <-
    failures %>%
    anti_join(names_resolved_to_other_sources, by = "query") %>%
    anti_join(names_with_mult_syns, by = "query")
  
  ### Check to make everything is accounted for
  assertthat::assert_that(isTRUE(all.equal(
    bind_rows(
      non_valid_queries %>% select(query),
      names_resolved_auto %>% select(query),
      names_with_mult_syns %>% select(query) %>% unique,
      names_resolved_to_other_sources %>% select(query),
      total_failures %>% select(query)
    ) %>% pull(query) %>% sort,
    name_resolution_results %>% pull(query) %>% sort
  )))
  
  list(
    names_resolved_auto = names_resolved_auto,
    names_resolved_to_other_sources = names_resolved_to_other_sources,
    names_with_mult_syns = names_with_mult_syns,
    total_failures = total_failures,
    non_valid_queries = non_valid_queries
  )
  
}

#' Take the results of automatic and manual name resolution and use
#' them to finalize GenBank data, subsetting to only resolved names
#'
#' @param names_resolved_auto Names that were successfully resolved automatically
#'   to Catalog of Life (i.e., World Ferns)
#' @param names_resolved_to_other_sources Names that were successfully resolved 
#'   automatically a source other than Catalog of Life
#' @param name_resolution_syns_to_use Names with multiple synonyms that have
#'   been manual inspected and the synonym to use selected
#' @param combined_metadata Tibble including sequences and associated
#' metadata. Sequences in column 'seq' as list-vector.
#'
#' @return Tibble with resolved names and sequences
#' 
resolve_genbank_names_final <- function (
  names_resolved_auto, names_resolved_to_other_sources, name_resolution_syns_to_use,
  combined_metadata) {
  
  # - names resolved automatically to CoL
  names_resolved_auto <- 
    names_resolved_auto %>%
    select(query, scientificName)
  
  # - names resolved automatically to other sources (GBIF, IPNI, Tropicos)
  names_resolved_to_other_sources <- 
    names_resolved_to_other_sources %>%
    select(query, scientificName = matched_name)
  
  # - names with mult. synonyms resolved manually
  name_resolution_syns_to_use <-
    name_resolution_syns_to_use %>% filter(!is.na(name_resolved_manual)) %>%
    select(query, scientificName = name_resolved_manual)
  
  names_resolved <- 
    bind_rows(
      names_resolved_auto, 
      names_resolved_to_other_sources, 
      name_resolution_syns_to_use) %>%
    assert(is_uniq, query) %>%
    taxastand::add_parsed_names(scientificName, species)
  
  # Update names in rbcL data to resolved names, and 
  # only keep sequences with resolved names.
  combined_metadata %>%
    rename(query = species) %>%
    inner_join(names_resolved, by = "query") %>%
    select(-query) %>%
    assert(not_na, scientificName, species)
  
}

#' Run all-by-all BLAST to detect rogue sequences in GenBank pteridophytes
#' 
#' All sequences will be BLASTed against each other.
#'
#' @param metadata_with_seqs Tibble containting sequences downloaded
#' from GenBank (column "seqs") and associated metadata (at least "accession"
#' and "species").
#' @param ppgi Tibble; genus-level and higher taxonomic system for pteridophytes,
#' according to PPGI 2016.
#' @param ... Additional arguments not used by this function but meant for 
#' tracking with drake
#' 
blast_rogues <- function (metadata_with_seqs, ...) {
  
  ### All-by-all BLAST ###
  
  # Create OTU column for naming sequences as species_accession
  metadata_with_seqs <- dplyr::mutate(
    metadata_with_seqs,
    otu = glue("{species}_{accession}") %>% stringr::str_replace_all(" ", "_")
  )
  
  # Extract sequences from metadata and rename
  pterido_seqs <- ape::as.DNAbin(metadata_with_seqs$seq)
  names(pterido_seqs) <- metadata_with_seqs$otu
  
  # Make sure that went OK
  assertthat::assert_that(is.list(pterido_seqs))
  assertthat::assert_that(inherits(pterido_seqs, "DNAbin"))
  assertthat::assert_that(all(names(pterido_seqs) == metadata_with_seqs$otu))
  
  # Blast to exclude rogues
  
  # Make temporary working dir for BLAST functions.
  blast_dir <- fs::dir_create(
    fs::path(
      tempdir(), 
      # Use hash as unique name for folder
      digest::digest(pterido_seqs)
    ))
  
  # Remove any gaps
  pterido_seqs <- ape::del.gaps(pterido_seqs)
  
  # Write out gapless sequences.
  ape::write.FASTA(pterido_seqs, fs::path(
    blast_dir, 
    "pterido_seqs.fasta"))
  
  # Make blast db.
  baitfindR::build_blast_db(
    in_seqs = fs::path(blast_dir, "pterido_seqs.fasta"),
    out_name = "pterido_seqs",
    wd = blast_dir,
    parse_seqids = FALSE
  )
  
  # Query all sequences
  baitfindR::blast_n(
    query = fs::path(blast_dir, "pterido_seqs.fasta"),
    database = "pterido_seqs",
    out_file = "blastn_results",
    other_args = c("-max_target_seqs", 10),
    wd = blast_dir,
    echo = TRUE
  )
  
  # Read in sequences.
  # BLAST doesn't output column headers, so we need to specify 
  # (make sure they match correctly first!)
  fmt6_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  
  blast_results <- readr::read_tsv(
    fs::path(blast_dir, "blastn_results"),
    col_names = fmt6_cols
  )
  
  # Cleanup
  fs::dir_delete (blast_dir)
  
  blast_results
  
}

#' Detect rogues from running all-by-all BLAST
#' 
#' Any sequence with a three top non-self hits to a different family 
#' will be flagged as a "rogue" sequence.
#' 
#' All families in Cyatheales and Saccolomatineae, respectively,
#' are treated as single families.
#'
#' @param metadata_with_seqs Tibble containting sequences downloaded
#' from GenBank (column "seqs") and associated metadata (at least "accession"
#' and "species").
#' @param blast_results Output of running blast_rogues()
#' @param ppgi Tibble; genus-level and higher taxonomic system for pteridophytes,
#' according to PPGI 2016.
#' @param ... Additional arguments not used by this function but meant for 
#' tracking with drake
#' 
#' @return Tibble; list of rogue sequences showing which families were matched
detect_rogues <- function(metadata_with_seqs, blast_results, ppgi, ...) {
  
  ### Detect rogues ###
  
  ### Check for monotypic families ###
  # These by default will match to a different (non-self)
  # family, so are not valid for checking rogues.
  exclude_from_rogues <-
    metadata_with_seqs %>% 
    dplyr::mutate(genus = stringr::str_split(species, " ") %>% purrr::map_chr(1)) %>%
    dplyr::left_join(
      select(ppgi, genus, family, suborder, order), 
      by = "genus") %>%
    dplyr::add_count(family) %>%
    dplyr::filter(n == 1) %>%
    dplyr::mutate(
      otu = glue("{species}_{accession}") %>% stringr::str_replace_all(" ", "_")
    ) %>%
    dplyr::pull(otu)
  
  # FIXME: need to run taxonomic resolution before rogue detection
  # e.g., NCBI name of Diplaziopsis cavaleriana (Christ) C. Chr. (Diplaziopsidaceae)
  # is Diplazium cavalerianum (Christ) M. Kato (Athyriaceae)
  #
  # Wait on this until we can use a fixed version of World Ferns though because
  # currenty it has many species placed in the wrong family, which will cause
  # more problems that it solves.
  
  # Group small (esp. monotypic) families
  # to avoid false-positives
  ppgi <- mutate(ppgi, family = case_when(
    order == "Cyatheales" ~ "Cyatheales",
    suborder == "Saccolomatineae" ~ "Saccolomatineae",
    TRUE ~ family
  ))
  
  # Make list of rogue sequences (accessions) to exclude
  # Start with blast results
  blast_results %>%
    # Remove extra '_R_' added by MAFFT when reversing seqs
    dplyr::mutate_at(vars(qseqid, sseqid), ~ stringr::str_remove_all(., "_R_")) %>%
    # Exclude monotypic families
    dplyr::filter(!qseqid %in% exclude_from_rogues) %>%
    # Filter out self-hits
    dplyr::filter(qseqid != sseqid) %>%
    # Keep top 3 best hits per query
    dplyr::group_by(qseqid) %>%
    dplyr::arrange(evalue, desc(bitscore), sseqid) %>%
    dplyr::slice(1:3) %>%
    dplyr::ungroup() %>%
    # Add columns for query genus and hit (i.e., 'subject') genus
    dplyr::mutate(
      q_genus = stringr::str_split(qseqid, "_") %>% purrr::map_chr(1),
      s_genus = stringr::str_split(sseqid, "_") %>% purrr::map_chr(1)
    ) %>%
    # Add columns for query family and hit family (look up from PPGI taxonomy)
    dplyr::left_join(
      select(ppgi, q_genus = genus, q_family = family),
      by = "q_genus"
    ) %>% 
    dplyr::left_join(
      select(ppgi, s_genus = genus, s_family = family),
      by = "s_genus"
    ) %>%
    # Consider rogue to be when top 3 hit families all the same and
    # different from the query family
    group_by(qseqid) %>%
    summarize(
      q_family = unique(q_family) %>% paste(collapse = ", "),
      s_family = unique(s_family) %>% paste(collapse = ", ")
    ) %>%
    filter(str_detect(s_family, ",", negate = TRUE)) %>%
    filter(q_family != s_family) %>%
    mutate(accession = str_split(qseqid, "_") %>% map_chr(last))
  
}

# Tracked version of anti-join for drake
anti_join_tracked <- function (id, ...) {
  
  dplyr::anti_join(...)
  
}

# Filter to one seq per species by removing 'rogue' sequences, then
# selecting the one
# with most "good" bases, then extract the sequences into a DNA object.
# Call ties by sorting on accession. Oviously, this
# shouldn't matter for the quality of sequence,
# but is reproducible.
filter_and_extract_pterido_rbcl <- function (metadata_with_seqs) {
  
  metadata_with_seqs <-
    metadata_with_seqs %>%
    # Only keep seqs > 99 bp
    filter(length > 99) %>%
    # Remove rogue sequences
    filter(is.na(rogue_summary)) %>%
    # Count number of missing nucleotides
    dplyr::mutate(
      n_missing = purrr::map_dbl(seq, ~ stringr::str_count(., "[^atgc]") %>% sum()),
      # The number of "good" bases is total length minus missing
      n_good_bases = slen - n_missing) %>%
    dplyr::group_by(species) %>%
    # Sort by most "good" bases
    dplyr::arrange(dplyr::desc(n_good_bases), accession) %>%
    # Take the top one
    dplyr::slice(1) %>%
    # Set OTU name for naming sequences
    dplyr::mutate(
      otu = glue("{species}_{accession}") %>% stringr::str_replace_all(" ", "_")
    )
  
  # Extract sequences from metadata and rename
  seqs <- ape::as.DNAbin(metadata_with_seqs$seq)
  names(seqs) <- metadata_with_seqs$otu
  
  # Make sure that went OK
  assertthat::assert_that(is.list(seqs))
  assertthat::assert_that(inherits(seqs, "DNAbin"))
  assertthat::assert_that(all(names(seqs) == metadata_with_seqs$otu))
  
  seqs
}

#' Fetch a specific gene from a list of accessions
#' 
#' Sometimes a given GenBank accession contains multiple genes or loci.
#' This grabs just the gene of interest from a given accession.
#' 
#' If the gene is not present in a given accession, that accession
#' will not appear in the results.
#'
#' @param accession Genbank accession number(s)
#' @param gene Target gene (must match gene name in "gene" field of
#' FEATURES exactly)
#'
#' @return List of class "DNAbin". The selected sequences.
#'
#' @examples 
#' fetch_gene_from_accessions("oi", "vey")
fetch_gene_from_accessions <- function (accession, gene) {
  
  assertthat::assert_that(assertthat::is.string(gene))
  
  assertthat::assert_that(is.character(accession))
  
  fetch_gene_from_genome_safely <- purrr::safely(gbfetch::fetch_gene_from_genome)
  
  first_try <- purrr::map(
    accession, 
    ~fetch_gene_from_genome_safely(gene = gene, accession = .)) %>%
    purrr::transpose
  
  success <- first_try$result %>% purrr::compact()
  
  failure <- first_try$error %>% purrr::compact()
  
  if(length(success) == 0) return (NULL)
  
  jntools::flatten_DNA_list(success) %>%
    purrr::set_names(names(.) %>% stringr::str_remove_all(paste0("-", gene)))
}

#' Tidy GenBank metadata
#'
#' GenBank metadata comes with a column "subtype" containing miscellaneous
#' data separated by '|'. The names of these data are in the "subname" column.
#' This function tidies these two columns, i.e., converts the data in the
#' "subtype" column so each value gets its own column.
#' 
#' @param data Dataframe; GenBank metadata obtained with
#' gbfetch::fetch_metadata.
#'
#' @return Dataframe.
#' 
#' @examples
#' raw_meta <- fetch_metadata("rbcl[Gene] AND Crepidomanes[ORGN]")
#' # Raw metadata still contains untidy data in the "subtype" and 
#' # "subname" columns
#' raw_meta
#' # Tidy!
#' tidy_genbank_metadata(raw_meta)

tidy_genbank_metadata <- function(data) {
  
  # Check assumptions: must have subtype and subname columns present
  assertthat::assert_that("subtype" %in% colnames(data))
  assertthat::assert_that("subname" %in% colnames(data))
  
  # Define helper function that works on one row at a time
  tidy_genbank_meta_single <- function(data) {
    
    # Early return if no data to parse in subtype
    if(data$subtype == "" | is.na(data$subtype)) return (data %>% dplyr::select(-subname, -subtype))
    if(data$subname == "" | is.na(data$subname)) return (data %>% dplyr::select(-subname, -subtype))
    
    # Split "subtype" into multiple columns
    # Use janitor::make_clean_names to de-duplicate names
    sub_cols <- data %>% 
      dplyr::pull(subtype) %>% 
      stringr::str_split("\\|") %>% 
      magrittr::extract2(1) %>% 
      janitor::make_clean_names()
    
    sub_data <- data %>% 
      dplyr::select(subname) %>% 
      tidyr::separate(subname, into = sub_cols, sep = "\\|")
    
    data %>% 
      dplyr::select(-subname, -subtype) %>% 
      dplyr::bind_cols(sub_data)
  }
  
  # Apply the function row-wise and combine back into a dataframe
  transpose(data) %>%
    purrr::map(as_tibble) %>%
    purrr::map_df(tidy_genbank_meta_single)
  
}

#' Filter GenBank sequences
#'
#' @param metadata_with_seqs Tibble containting sequences downloaded
#' from GenBank (column "seqs") and associated metadata (at least "accession"
#' and "species").
#' @param min_len Minimum sequence length to retain
#' @param ppgi Tibble; genus-level and higher taxonomic system for pteridophytes,
#' according to PPGI 2016.
#' @param ... Additional arguments not used by this function but meant for 
#' tracking with drake
#'
#' @return Tibble
filter_genbank_seqs <- function (metadata_with_seqs, min_len, ppgi, ...) {
  
  # Check for genera missing from ppgi, warn before dropping
  genus_not_in_ppgi <-
    metadata_with_seqs %>%
    dplyr::mutate(genus = stringr::str_split(species, " ") %>% purrr::map_chr(1)) %>%
    dplyr::anti_join(ppgi, by = "genus") %>%
    dplyr::mutate(msg = glue("{accession} ({genus})")) %>%
    dplyr::arrange(genus)
  
  assertthat::validate_that(
    nrow(genus_not_in_ppgi) == 0,
    msg = warning(glue(
      "{nrow(genus_not_in_ppgi)} accessions do not have genera in ppgi and will be dropped: 
       {paste(genus_not_in_ppgi$msg, collapse = ', ')}"))
  )
  
  # Apply filters
  metadata_with_seqs %>%
    dplyr::anti_join(genus_not_in_ppgi, by = "accession") %>%
    dplyr::filter(length > min_len) 
  
}

#' Select a set of GenBank genes by joining on common vouchers
#' 
#' This assumes genes comprise: "rbcL", "atpA", "atpB", "rps4"
#' 
#' @param genbank_seqs_tibble Tibble; metadata for GenBank sequences
#' including columns for `gene`, `species`, `specimen_voucher`, 
#' `accession` (GenBank accession), and `length` (gene length in bp)
#' @param n_seqs_per_sp String; format of output data. If 'single', only one set 
#' of sequences per species will be returned, prioritizing species with rbcL. 
#' If 'multiple', all sequences will be returned (multiple seqs per species)
#'
#' @return Tibble in wide format joining genes based on species + voucher
#' 
select_genbank_genes <- function (genbank_seqs_tibble, n_seqs_per_sp = c("single", "multiple")) {
  
  assertthat::assert_that(assertthat::is.string(n_seqs_per_sp))
  
  assertthat::assert_that(
    n_seqs_per_sp %in% c("single", "multiple"),
    msg = "'n_seqs_per_sp' must be either 'single' or 'multiple'")
  
  ### Join genes based on voucher ###
  # Original genbank data is in "long" format, with one row per accession
  # We need to convert this to wide by joining gene sequences with the same voucher.
  # To join without combining NA vouchers, split into a list of dataframes,
  # then do a full join with na_matches = "never"
  gb_data <- 
    genbank_seqs_tibble %>%
    select(gene, species, specimen_voucher, accession, length) %>%
    # Select one sequence per specimen per species per gene based on gene length
    group_by(gene, specimen_voucher, species) %>%
    arrange(desc(length)) %>%
    slice(1) %>%
    ungroup %>%
    # FIXME: work-around to avoid joining on NA vouchers: treat each as a distinct sample
    mutate(
      row_num = 1:nrow(.),
      specimen_voucher = ifelse(is.na(specimen_voucher), glue("specimen_missing_{row_num}"), specimen_voucher)
    ) %>%
    # Convert to wide format, joining on voucher
    # - first split into a list of dataframes by gene
    group_by(gene) %>%
    group_split %>%
    # For each dataframe, convert to wide format and
    # rename "length" and "accession" columns by gene
    map(
      ~pivot_wider(.,
                   id_cols = c("species", "specimen_voucher"),
                   names_from = gene,
                   values_from = c("length", "accession")
      )
    ) %>%
    # Join the gene sequences by species + voucher
    # FIXME: na_matches = "never" is currently breaking with full_join()
    # once this bug gets fixed in dplyr, can delete the work-around above
    reduce(full_join, by = c("species", "specimen_voucher"), na_matches = "na") %>%
    mutate_at(vars(contains("length")), ~replace_na(., 0)) %>%
    # Add total length of all genes (need to sum row-wise)
    # https://stackoverflow.com/questions/31193101/how-to-do-rowwise-summation-over-selected-columns-using-column-index-with-dplyr
    mutate(total_length = pmap_dbl(select(., contains("length")), sum)) %>%
    # FIXME: last part of work-around: convert "specimen_missing" back to NA. remove this once bug gets fixed in dplyr.
    mutate(specimen_voucher = ifelse(str_detect(specimen_voucher, "specimen_missing"), NA, specimen_voucher))
  
  # If `n_seqs_per_sp` is 'multiple' exit early, returning all seqs (including multiple seqs per species) 
  if (n_seqs_per_sp == "multiple") return (gb_data)
  
  # Otherwise continue on, selecting one set of sequences per species.
  
  ### Select final sequences ###
  
  # Highest priority species: those with rbcL and at least one other gene.
  # Choose best specimen per species by total sequence length
  rbcl_and_at_least_one_other <-
    gb_data %>% 
    filter(!is.na(accession_rbcL)) %>%
    filter_at(
      vars(accession_atpA, accession_atpB, accession_rps4), 
      any_vars(!is.na(.))) %>%
    group_by(species) %>%
    arrange(desc(total_length)) %>%
    slice(1) %>%
    ungroup
  
  # Next priority: anything with rbcL only
  rbcl_only <-
    gb_data %>%
    filter(!species %in% rbcl_and_at_least_one_other$species) %>%
    filter(!is.na(accession_rbcL)) %>%
    group_by(species) %>%
    arrange(desc(total_length)) %>%
    slice(1) %>%
    ungroup
  
  # Next priority: any other species on basis of total seq. length
  other_genes <-
    gb_data %>%
    filter(!species %in% rbcl_and_at_least_one_other$species) %>%
    filter(!species %in% rbcl_only$species) %>%
    group_by(species) %>%
    arrange(desc(total_length)) %>%
    slice(1) %>%
    ungroup
  
  # Combine into final list
  bind_rows(
    rbcl_and_at_least_one_other,
    rbcl_only,
    other_genes
  )
  
}

# Helper function to extract DNA sequences from tibble
# where they are stored as a list-column
extract_seqs <- function(data) {
  seqs <- data$seq
  names(seqs) <- data$accession
  ape::as.DNAbin(seqs)
}

#' Rename genes from number (1, 2, etc) to gene name ("rbcL", "atpA", etc)
#'
#' @param raw_fasta_all_genes Combined dataframe including metadata
#' and DNA sequence for all downloaded GenBank (Sanger) sequences
#' @param genes_used Character vector of genes used
#'
#' @return Dataframe with genes listed by name instead of number
#' 
rename_genes <- function (raw_fasta_all_genes, genes_used) {
  
  ### Recode gene names from numbers to actual genes used ###
  # Note this depends on the order of `genes_used`:
  # must be in the same order as `target_used used` for fetch_fern_gene()
  genes_used <- genes_used %>% set_names(1:length(genes_used))
  
  raw_fasta_all_genes %>%
    rename(gene_number = gene) %>%
    mutate(
      gene = recode(gene_number, !!!genes_used)
    )
  
}

#' Filter out species in plastome data from Sanger data
#'
#' @param plastid_genes_unaligned List of unaligned plastome genes
#' @param plastome_metadata_renamed Plastome metadata after standardizing names
#' @param genbank_accessions_selection Selected GenBank (Sanger) sequences to use
#' @param filter Logical; should filtering be done or not?
#'
#' @return genbank_accessions_selection, with species in plastome data filtered out (if `filter` = TRUE) 
#' 
filter_out_plastome_species <- function (plastid_genes_unaligned, plastome_metadata_renamed, genbank_accessions_selection, filter) {
  
  # Don't do any filtering if `filter` is false
  if(filter == FALSE) return (genbank_accessions_selection)
  
  ### Make list of species in selected plastome genes ###
  # (not the same as plastome_metadata_renamed, since that was ALL plastomes, and
  # we need a list of species names of only the plastomes that will actually be used)
  plastid_genes_unaligned_species_list <-
    map_df(plastid_genes_unaligned, ~names(.) %>% tibble(accession = .), .id = "gene") %>%
    left_join(plastome_metadata_renamed, by = "accession") %>%
    select(gene, accession, species) %>%
    assert(not_na, species)
  
  # Filter out species from GenBank (Sanger) accessions that are already in plastome data
  anti_join(
    genbank_accessions_selection, 
    plastid_genes_unaligned_species_list,
    by = "species"
  )
  
}

#' Combine sequences for target genes from GenBank with sequences
#' extracted from plastomes
#'
#' @param raw_fasta_all_genes Tibble of DNA sequence data downloaded from GenBank
#' @param genbank_accessions_selection Final selection of accessions to use after removing rogues
#' @param plastid_genes_unaligned List of unaligned genes extracted from plastomes
#'
#' @return List of class DNAbin
#' 
combine_genbank_with_plastome <- function (
  raw_fasta_all_genes,
  genbank_accessions_selection,
  plastid_genes_unaligned
) {
  
  # Convert final selected GenBank (Sanger) accessions to long format
  final_gb_accessions <-
    genbank_accessions_selection %>%
    select(species, specimen_voucher, contains("accession")) %>%
    pivot_longer(
      cols = contains("accession"),
      names_to = "gene",
      names_pattern = "_(.*)$",
      values_to = "accession") %>%
    filter(!is.na(accession))
  
  # GenBank (Sanger) data: add sequence data, group by gene
  final_seqs_grouped <-
    final_gb_accessions %>%
    # Check that the combination of gene + accession is unique
    assert_rows(col_concat, is_uniq, accession, gene) %>%
    left_join(
      select(raw_fasta_all_genes, gene, accession, seq),
      # Join on gene + accession, since some different genes share the same acc
      by = c("gene", "accession")
    ) %>%
    # Check that the combination of gene + accession is unique
    assert_rows(col_concat, is_uniq, accession, gene) %>%
    # Set grouping
    group_by(gene)
  
  # GenBank (Sanger) data: convert to list of DNA sequences, name by gene
  genbank_genes_unaligned <-
    final_seqs_grouped %>%
    group_split %>%
    map(extract_seqs)
  
  names(genbank_genes_unaligned) <- group_keys(final_seqs_grouped) %>% pull(gene)
  
  ### Combine with plastome sequences ###
  
  # - Make vector of genes in common between Sanger and plastome sequences
  common_gene_names <- intersect(names(genbank_genes_unaligned), names(plastid_genes_unaligned)) 
  
  # - Make a list of accessions for genes in common between Sanger and plastome sequences
  common_genes <-
    common_gene_names %>%
    map(~c(genbank_genes_unaligned[[.]], plastid_genes_unaligned[[.]])) %>%
    set_names(common_gene_names)
  
  # - Make a list of accessions in Sanger genes only (likely zero, but for completeness' sake)
  genbank_genes_only <- genbank_genes_unaligned %>%
    magrittr::extract(setdiff(names(genbank_genes_unaligned), names(plastid_genes_unaligned)))
  
  # - Make a list of accessions in plastome sequences only
  plastome_genes_only <- plastid_genes_unaligned %>%
    magrittr::extract(setdiff(names(plastid_genes_unaligned), names(genbank_genes_unaligned)))
  
  # Combine the lists
  c(common_genes, genbank_genes_only, plastome_genes_only)
  
}

#' Make a voucher lookup table
#'
#' @param genbank_seqs_names_resolved Metadata of plastid (Sanger) sequences
#' from Genbank, with names resolved
#' @param genbank_accessions_selection_multiple Selected plastid (Sanger) accessions
#' to use for tree with multiple tips per species
#' @param plastid_selection_voucher Selected plastome accessions
#' to use for tree with multiple tips per species
#'
#' @return Dataframe
#' 
make_voucher_table <- function(
  genbank_seqs_names_resolved,
  genbank_accessions_selection_multiple,
  plastid_selection_voucher
) {
  
  # Extract list of all accessions used in the final genbank accessions selection
  # (multiple sequences per species)
  # - Sanger sequences: this has already been filtered to one sequence per specimen per gene
  all_accessions_sanger <- genbank_accessions_selection_multiple %>%
    select(contains("accession")) %>%
    unlist %>%
    magrittr::extract(!is.na(.))
  
  # Plastome sequences: these are also already selected (best plastome per voucher per species)
  all_accessions_plastome <- plastid_selection_voucher$accession %>% unique
  
  all_accessions <- c(all_accessions_sanger, all_accessions_plastome) %>% unique
  
  # Double check that all vouchers appear no more than once per gene
  genbank_seqs_names_resolved %>%
    select(gene, accession, species, voucher = specimen_voucher) %>% 
    unique %>%
    filter(accession %in% all_accessions) %>%
    group_by(species, voucher) %>%
    count(voucher, sort = TRUE) %>%
    ungroup %>%
    assert(within_bounds(1, length(target_genes)), n, error_fun = assertr::error_stop, success_fun = assertr::success_logical)
  
  # Make voucher-lookup table: assign a voucher ID to each accession within each species/gene combination
  
  # First, combine plastome and sanger data
  bind_rows(
    plastid_selection_voucher %>% 
      select(gene, accession, species, voucher = specimen_voucher),
    genbank_seqs_names_resolved %>%
      select(gene, accession, species, voucher = specimen_voucher)
  ) %>% 
    # Add "voucher ID" number for each voucher within a species
    unique %>%
    # - if no voucher available, consider each accession unique
    mutate(voucher = ifelse(is.na(voucher), accession, voucher)) %>%
    filter(accession %in% all_accessions) %>%
    group_by(species) %>%
    mutate(voucher_id = factor(voucher) %>% as.numeric) %>%
    ungroup %>%
    # Make sure no missing values (except for voucher)
    assert(not_na, gene, accession, species, voucher_id) %>%
    # Make sure combination of gene + species + voucher_id is unique
    assert_rows(col_concat, is_uniq, gene, species, voucher_id, error_fun = assertr::error_stop) %>%
    # Rename species as species plus voucher ID
    mutate(species = paste(species, voucher_id))
  
}

# Download plastid genes from plastomes ----

#' Download plastome metadata
#'
#' @param start_date Earliest date to download
#' @param end_date Most recent date to download
#' @param outgroups Dataframe of genbank accession numbers
#' and species names to use for outgroups (taxa that wouldn't
#' be captured by the query)
#'
#' @return Tibble of GenBank metadata combining the results of
#' the query and the outgroups.
#' 
download_plastome_metadata <- function (start_date = "1980/01/01", end_date, outgroups) {
  
  assertthat::assert_that(assertthat::is.string(start_date))
  
  assertthat::assert_that(assertthat::is.string(end_date))
  
  # Format GenBank query: all ferns plastomes within specified dates
  query = glue('genome AND Polypodiopsida[ORGN] AND (plastid OR chloroplast) AND (partial OR complete) AND ("{start_date}"[PDAT]:"{end_date}"[PDAT])')
  
  # Download metadata and accession numbers.
  pterido_metadata <- gbfetch::fetch_metadata(query) %>%
    # "genome" hit in some unexpeted places, so filter out short seqs
    # that are not actually plastomes
    filter(slen > 7000) %>%
    # Filter out a strangely formatted accession that causes errors
    # with biofiles::getFeatures() and biofiles::filter()
    # (AP004638, Psilotum nudum)
    filter(accession != "AP004638") %>%
    # Parse out data contained in the "subtype" and "subname" columns
    tidy_genbank_metadata() %>%
    # Remove weird brackets from species names
    mutate(species = str_remove_all(species, "\\[") %>% str_remove_all("\\]")) %>%
    # Remove GenBank duplicate sequences starting with "NC_"
    # - consider accessions with same species and seq length to be same
    # if they only differ in accession.
    mutate(maybe_dup = case_when(
      str_detect(accession, "NC_") ~ 1,
      TRUE ~ 0)) %>%
    group_by(species, slen) %>%
    arrange(maybe_dup) %>%
    slice(1) %>%
    ungroup %>%
    # Sort by accession
    arrange(accession)
  
  # Download outgroup metadata: use one representative of each major 
  # seed plant, lycophyte, and bryo group
  og_query <- outgroups %>% pull(accession) %>% paste(collapse = "[accession] OR ") %>%
    paste("[accession]", collapse = "", sep = "")
  
  outgroup_metadata <- gbfetch::fetch_metadata(og_query) %>%
    tidy_genbank_metadata() 
  
  bind_rows(pterido_metadata, outgroup_metadata) %>%
    select(-maybe_dup)
  
}

#' Fetch a set of genes from a plastome
#'
#' @param genes Vector of gene names
#' @param accession Genbank accession number (of a plastome)
#' @param max_length Maximum length to accept for genes. Used to filter
#' out any abnormally (probably erroneously) long genes.
#'
#' @return List of lists; each list is a gene sequence of class DNAbin
#' with a single sequence.
#' The lists are named by gene, and the sequences are named after the accession.
#' 
#' test <- fetch_fern_genes_from_plastome(
#' accession = "AY178864",
#' max_length = 10000,
#' genes = read_lines("data_raw/wei_2017_coding_genes.txt")
#' )
fetch_fern_genes_from_plastome <- function (genes, accession, max_length = 10000) {
  
  # Get GenBank ID for the accession
  uid <- reutils::esearch(term = accession, db = "nucleotide", usehistory = TRUE)
  
  # Make sure there is only 1 hit for that accession
  num_hits <- reutils::content(uid, as = "text") %>% str_match("<eSearchResult><Count>([:digit:]+)<\\/Count>") %>% magrittr::extract(,2)
  
  assertthat::assert_that(
    num_hits == 1,
    msg = "Did not find exactly one accession")
  
  # Download complete GenBank record and write it to a temporary file
  temp_dir <- tempdir()
  temp_file <- fs::path(temp_dir, "gb_records.txt")
  
  reutils::efetch(uid, "nucleotide", rettype = "gb", retmode = "text", outfile = temp_file)
  
  # Parse flatfile
  gb_entry <- readr::read_file(temp_file)
  
  # Some genes are missing or duplicates for a given plastome,
  # so avoid errors by wrapping extract_sequence() in safely()
  extract_sequence_safely <- safely(extract_sequence)
  
  # get the results
  extracted_genes <- map(genes, ~extract_sequence_safely(gb_entry, .)) %>%
    transpose() %>%
    pluck("result")
  
  # subset gene names to those to successfully extracted
  genes_successful <- genes[!map_lgl(extracted_genes, is.null)]
  
  # Subset results to successful genes, filter by length, and set names
  extracted_genes_filtered <-
    extracted_genes %>%
    # Drop errors
    compact() %>%
    flatten() %>%
    # Name each list item as the gene
    set_names(genes_successful) %>%
    # Convert to DNAbin
    map(~stringr::str_split(., "") %>% ape::as.DNAbin()) %>%
    # Exclude abnormally long sequences
    # (rps12 from Adiantum capillus-veneris AY178864) with 72,969 bp!
    # Appears to be an error in annotation?
    magrittr::extract(!map_lgl(., ~map_dbl(., length) > max_length)) %>%
    # Name each DNAbin as the accession
    map(~set_names(., accession))
  
  # Return as a list of 1, so the results can be combined into a list later
  list(extracted_genes_filtered) %>% set_names(accession)
  
}

#' Helper function to filter accessions missing > 50% of genes
#' and genes absent from > 50% of sequences used in select_plastid_seqs()
#'
#' @param seq_list List of sequences
#'
#' @return Length of each sequence in long format
get_gene_lengths <- function (seq_list) {
  seq_list %>%
    map_df(~map_dbl(., length), .id = "gene") %>%
    rename(slen = 2)
}

#' Helper function to to filter accessions missing > 50% of genes
#' and genes absent from > 50% of sequences used in select_plastid_seqs()
#'
#' @param gene_lengths_best Dataframe of 'best' (ie, least missing data)
#' plastid genes
#'
#' @return Dataframe
#'
filter_majority_missing <- function (gene_lengths_best) {
  
  # Helper function to to filter accessions missing > 50% of genes
  # and genes absent from > 50% of sequences.
  
  # Assume gene is missing if number of base pairs is 2 or less
  # (not sure why but there are few of these, then most genes have at least ~40 bp)
  missing_acc_summ <-
    gene_lengths_best %>%
    mutate(missing = case_when(
      slen < 2 ~ TRUE,
      TRUE ~ FALSE
    )) %>%
    group_by(gene) %>%
    summarize(
      n_missing_accs = sum(missing),
      .groups = "drop"
    ) %>%
    mutate(
      keep_gene = case_when(
        n_missing_accs < 0.5 * n_distinct(gene_lengths_best$accession) ~ TRUE,
        TRUE ~ FALSE
      )
    )
  
  missing_gene_summ <-
    gene_lengths_best %>%
    mutate(missing = case_when(
      slen < 2 ~ TRUE,
      TRUE ~ FALSE
    )) %>%
    group_by(accession) %>%
    summarize(
      n_missing_genes = sum(missing),
      .groups = "drop"
    ) %>%
    mutate(
      keep_acc = case_when(
        n_missing_genes < 0.5 * n_distinct(gene_lengths_best$gene) ~ TRUE,
        TRUE ~ FALSE
      )
    )
  
  gene_lengths_best %>%
    left_join(missing_acc_summ, by = "gene") %>%
    left_join(missing_gene_summ, by = "accession") %>%
    filter(keep_gene, keep_acc)
}

#' Select a list of plastid sequences to use from a list of plastid genes and
#' plastome metadata
#'
#' @param plastid_seq_list List of plastid genes. Each list item is
#' a named list of gene sequences for a plastome accession.
#' @param plastome_metadata Associated plastome metadata (species names)
#' @param filter_by Should the list be selected by the best representative
#' accession per species, genus, or voucher?
#'
select_plastid_seqs <- function (plastid_seq_list, plastome_metadata, filter_by = c("species", "genus", "voucher")) {
  
  assertthat::assert_that(assertthat::is.string(filter_by))
  
  assertthat::assert_that(filter_by %in% c("species", "genus", "voucher"))
  
  # Drop any plastome sequences with zero genes
  check_genes <- function(plastome_seq) {
    length(plastome_seq) > 0
  }
  
  plastid_seq_list <- plastid_seq_list[map_lgl(plastid_seq_list, check_genes)]
  
  # Make tibble of gene lengths by accession, including species and voucher
  gene_lengths <- 
    plastid_seq_list %>%
    map_df(get_gene_lengths, .id = "plastid_seq_name") %>%
    mutate(accession = str_remove(plastid_seq_name, "clean_plastid_seqs_plastid_seqs_")) %>%
    left_join(select(plastome_metadata, accession, species, specimen_voucher), by = "accession") %>%
    assert(not_na, species, gene, accession, slen)
  
  # Missing genes (length 0) are not in the original sequences list,
  # so add these by crossing all combinations of accession and gene
  gene_lengths <- 
    purrr::cross_df(list(
      gene = gene_lengths$gene %>% unique, 
      accession = gene_lengths$accession %>% unique)) %>%
    left_join(select(gene_lengths, gene, accession, slen), by = c("gene", "accession")) %>%
    left_join(select(gene_lengths, accession, species) %>% unique, by = "accession") %>%
    left_join(select(gene_lengths, accession, specimen_voucher) %>% unique, by = "accession") %>%
    left_join(select(gene_lengths, accession, plastid_seq_name) %>% unique, by = "accession") %>%
    mutate(slen = replace_na(slen, 0)) %>%
    assert(not_na, gene, accession, slen, species)
  
  # Get table of maximum lengths per gene
  # (we will assume these are the actual max. lengths)
  max_lengths <-
    gene_lengths%>%
    arrange(desc(slen)) %>%
    group_by(gene) %>%
    summarize(
      max_length = max(slen),
      .groups = "drop"
    )
  
  # Identify the "best" accessions as those having the least
  # amount of missing data overall per species
  best_accessions_by_species <-
    gene_lengths %>%
    left_join(max_lengths, by = "gene") %>%
    mutate(rel_len = slen / max_length) %>%
    assert(not_na, accession, species) %>%
    # first get total length for each accession
    group_by(accession, species) %>%
    summarize(
      total_rel_len = sum(rel_len),
      .groups = "drop"
    ) %>%
    # then sort by species and keep the one with the greatest length
    group_by(species) %>%
    arrange(desc(total_rel_len)) %>%
    slice(1) %>%
    ungroup
  
  # Make a table of (relative) gene lengths for
  # the best accession per species
  gene_lengths_best_by_species <-
    best_accessions_by_species %>%
    select(accession) %>%
    left_join(gene_lengths, by = "accession") %>%
    left_join(max_lengths, by = "gene") %>%
    mutate(rel_len = slen / max_length)
  
  # Do the same at the genus level:
  # Best accession per genus
  best_accessions_by_genus <-
    gene_lengths %>%
    left_join(max_lengths, by = "gene") %>%
    mutate(rel_len = slen / max_length) %>%
    # add Genus column
    mutate(genus = str_split(species, " ") %>% map_chr(1)) %>%
    # first get total length for each accession
    group_by(accession, genus) %>%
    summarize(
      total_rel_len = sum(rel_len),
      .groups = "drop"
    ) %>%
    ungroup %>%
    # then sort by genus and keep the one with the greatest length
    group_by(genus) %>%
    arrange(desc(total_rel_len)) %>%
    slice(1) %>%
    ungroup
  
  # And best gene lengths by genus
  gene_lengths_best_by_genus <-
    best_accessions_by_genus %>%
    select(accession) %>%
    left_join(gene_lengths, by = "accession") %>%
    left_join(max_lengths, by = "gene") %>%
    mutate(rel_len = slen / max_length)
  
  # Do the same at the voucher level:
  # Best accession per voucher
  best_accessions_by_voucher <-
    gene_lengths %>%
    left_join(max_lengths, by = "gene") %>%
    mutate(rel_len = slen / max_length) %>%
    # If voucher is NA (often the case), use accession as stand-in
    # (so assume each accession comes from a different voucher)
    mutate(specimen_voucher = ifelse(is.na(specimen_voucher), accession, specimen_voucher)) %>%
    assert(not_na, accession, specimen_voucher) %>%
    group_by(accession, specimen_voucher) %>%
    summarize(
      total_rel_len = sum(rel_len),
      .groups = "drop"
    ) %>%
    group_by(specimen_voucher) %>%
    arrange(desc(total_rel_len)) %>%
    slice(1) %>%
    ungroup
  
  # Make a table of (relative) gene lengths for
  # the best accession per voucher
  gene_lengths_best_by_voucher <-
    best_accessions_by_voucher %>%
    select(accession) %>%
    left_join(gene_lengths, by = "accession") %>%
    left_join(max_lengths, by = "gene") %>%
    mutate(rel_len = slen / max_length)
  
  # Select gene lengths filtered by genus, species, or voucher
  gene_lengths_best_filtered_by_genus <-
    filter_majority_missing(gene_lengths_best_by_genus) %>%
    select(accession, species, gene) %>%
    left_join(select(gene_lengths, plastid_seq_name, accession) %>% unique, by = "accession") %>%
    assert(not_na, everything())
  
  gene_lengths_best_filtered_by_species <-
    filter_majority_missing(gene_lengths_best_by_species)  %>%
    select(accession, species, gene) %>%
    left_join(select(gene_lengths, plastid_seq_name, accession) %>% unique, by = "accession") %>%
    assert(not_na, everything())
  
  gene_lengths_best_filtered_by_voucher <-
    filter_majority_missing(gene_lengths_best_by_voucher)  %>%
    select(accession, species, specimen_voucher, gene) %>%
    left_join(select(gene_lengths, plastid_seq_name, accession) %>% unique, by = "accession") %>%
    assert(not_na, accession, species, gene, plastid_seq_name)
  
  switch(
    filter_by,
    genus = gene_lengths_best_filtered_by_genus,
    species = gene_lengths_best_filtered_by_species,
    voucher = gene_lengths_best_filtered_by_voucher)
  
}

#' Reformat a list of plastid genes from 
#' a named list of gene sequences by plastome accession to
#' a named list of sequences from various accessions by gene
#'
#' @param plastid_seq_list  List of plastid genes. Each list item is
#' a named list of gene sequences for a plastome accession.
#' @param plastid_selection Tibble of selected plastomes to use.
#'
#' @return List of unaligned genes.
#' 
extract_seqs_by_gene <- function (plastid_seq_list, plastid_selection) {
  plastid_seq_list %>%
    magrittr::extract(unique(plastid_selection$plastid_seq_name)) %>%
    transpose %>%
    map(~purrr::compact(.) %>% purrr::reduce(c))
}

#' Run trimal on a DNA alignment in "auto" mode
#' 
#' trimal removes low-quality (i.e., poorly aligned, gappy sites) from
#' an alignment.
#'
#' @param seqs DNA sequence alignment; matrix of class 'DNAbin'
#' @param echo Optional; should the output of trimal be printed to the
#' screen?
#' @param ... Extra arguments not used by this function but meant for
#' tracking with drake.
#'
#' @return DNA sequence alignment; matrix of class 'DNAbin' with gappy sites
#' removed by trimal
#' 
trimal_auto <- function (seqs, echo = FALSE, ...) {
  
  # Check arguments
  assertthat::assert_that(inherits(seqs, "DNAbin"), msg = "seqs must be of class DNAbin")
  assertthat::assert_that(is.logical(echo))
  
  if(!is.matrix(seqs)) seqs <- as.matrix(seqs)
  
  # Check that trimal is installed and on the PATH
  check_trimal <- tryCatch({
    processx::run("trimal", "-h", echo = FALSE)
  }, warning = function(w) {
    stop("trimal not installed and on path")
  }, error = function(e) {
    stop("trimal not installed and on path")
  }, finally = {
    TRUE
  })
  
  # Write out alignment to temp dir
  temp_wd <- tempdir()
  
  in_file_name <- glue("{digest::digest(seqs)}.fasta")
  
  out_file_name <- glue("{digest::digest(seqs)}.aln")
  
  ape::write.FASTA(seqs, fs::path(temp_wd, in_file_name))
  
  # Set up arguments to trimal
  args <- c("-in", in_file_name, "-out", out_file_name, "-fasta", "-automated1")
  
  # Run trimal
  results <- processx::run("trimal", args, wd = temp_wd, 
                           echo = echo)
  
  # Read in results
  results <- ape::read.FASTA(fs::path(temp_wd, out_file_name)) %>% as.matrix
  
  # Clean up
  fs::file_delete(fs::path(temp_wd, out_file_name))
  fs::file_delete(fs::path(temp_wd, in_file_name))
  
  # All done
  results
}

#' Concatenate a list of aligned genes
#'
#' @param dna_list List of matrices of class DNAbin
#'
#' @return Matrix of class DNAbin
#'
#' @examples
#' data(woodmouse)
#' gene_1 <- woodmouse[,1:100]
#' gene_2 <- woodmouse[,101:200]
#' woodmouse_genes <- list(gene_1, gene_2)
#' concatenate_genes(woodmouse_genes)
concatenate_genes <- function (dna_list) {
  require(apex)
  
  assertthat::assert_that(is.list(dna_list))
  assertthat::assert_that(all(lapply(dna_list, class) == "DNAbin"), 
                          msg = "All elements of dna_list must be of class DNAbin")
  assertthat::assert_that(all(sapply(dna_list, is.matrix)), 
                          msg = "All elements of dna_list must be matrices")
  
  # Check that there are no duplicate sequence names (species) within a gene
    map_df(dna_list, ~rownames(.) %>% tibble(species = .), .id = "gene") %>%
      assert_rows(col_concat, is_uniq, species, gene, error_fun = assertr::error_stop)
  
  dna_multi <- new("multidna", dna_list) 
  apex::concatenate(dna_multi)
}

#' Resolve species names in plastome metadata using Catalog of Life
#' as the taxonomic standard
#'
#' @param plastome_metadata Dataframe with column "species" containing
#' original species names from GenBank
#' @param col_plants Dataframe; taxonomic standard based on the
#' Catalog of Life
#'
#' @return Dataframe; plastome_metadata with new column `resolved_name`
#' containing the standardized name attached to it.
#' 
resolve_pterido_plastome_names <- function(plastome_metadata, col_plants) {
  
  # Resolve names in plastome data:
  # - first automatically
  name_resolution_results <- taxastand::resolve_fern_names(
    names = unique(plastome_metadata$species),
    col_plants = col_plants,
    exact_match_by = "species",
    resolve_to = "scientific_name")
  
  # - Manual fixes for any names that couldn't be resolved manually
  name_resolution_results_fixed <-
    name_resolution_results %>%
    mutate(scientificName = case_when(
      exclude_non_pterido_genus == TRUE ~ query,
      query == "Physcomitrella patens" ~ query,
      query == "Physcomitrium patens" ~ query,
      query == "Marchantia polymorpha" ~ query,
      query == "Psilotum nudum" ~ query,
      query == "Botrychium sp. ternatum/japonicum" ~ "Sceptridium ternatum",
      query == "Adiantum shastense" ~ query,
      query == "Cryptogramma acrostichoides" ~ query,
      query == "Pecluma dulcis" ~ query,
      query == "Alsophila podophylla" ~ query,
      query == "Asplenium nidus" ~ query,
      query == "Athyrium sinense" ~ query,
      query == "Dicksonia antarctica" ~ query,
      query == "Diplazium striatum" ~ query,
      query == "Lygodium microphyllum" ~ query,
      TRUE ~ scientificName)) %>%
    # Make sure all names are accounted for
    assert(not_na, scientificName) %>%
    # Add just species name
    taxastand::add_parsed_names(scientificName, species) %>%
    # Manual fix after gnparser drops "nudum" from Psilotum nudum
    mutate(species = case_when(
      species == "Psilotum" ~ "Psilotum nudum",
      TRUE ~ species
    )) %>%
    # Make sure all species have genus and epithet separated by space
    assert(function (x) str_detect(x, " "), species) %>%
    select(query, resolved_name = species)
  
  plastome_metadata %>%
    left_join(select(name_resolution_results_fixed, species = query, resolved_name), by = "species") %>%
    mutate(species = resolved_name) %>%
    assert(not_na, species) %>%
    select(-resolved_name)
}

#' Combine GenBank rbcL and plastome-derived rbcL sequences
#' 
#' Preferentially use plastome-derived rbcL sequences if the same
#' species is present in both.
#' 
#' All spaces in species names replaced with underscores (for proper
#' handling in phy. analysis)
#'
#' @param plastid_genes_unaligned Unaligned plastid genes from plastomes
#' @param plastome_metadata_renamed Plastome metadata with resolved names
#' @param pterido_rbcl_clean_seqs Unaligned rbcL sequences from genbank with
#' resolved names
#'
#' @return List; first item is the DNA sequences (list of class DNAbin); second
#' item is a dataframe with GenBank accession numbers for each species.
#' @export
#'
#' @examples
combine_rbcL <- function(plastid_genes_unaligned, plastome_metadata_renamed, pterido_rbcl_clean_seqs) {
  
  # Extract rbcL from the plastome genes list, and rename by species
  plastome_rbcL <- plastid_genes_unaligned[["rbcL"]]
  
  plastome_rbcL_names <-
    tibble(accession = names(plastome_rbcL)) %>%
    left_join(
      select(plastome_metadata_renamed, accession, species),
      by = "accession") %>%
    mutate(species = str_replace_all(species, " ", "_"))
  
  names(plastome_rbcL) <- plastome_rbcL_names$species
  
  # Make tibble of rbcL alignment tip names, species, accessions,
  # and column if that species is in the plastome data
  rbcL_species <-
    tibble(tip_name = names(pterido_rbcl_clean_seqs)) %>%
    # Make sure there aren't any accessions with an underscore like 'NC_178345'
    # because these would break the str_match used to extract accession numbers
    assert(function (x) str_detect(x, "[:upper:]_", negate = TRUE), tip_name) %>%
    mutate(
      species = map_chr(tip_name, ~str_match(., "^(.*)_.*$") %>% magrittr::extract(,2)),
      accession = map_chr(tip_name, ~str_match(., "^.*_(.*)$") %>% magrittr::extract(,2)),
      in_plastomes = species %in% names(plastome_rbcL)
    )
  
  # Make list of rbcL seqs to exclude because those species are in
  # the plastome data
  rbcL_to_exclude <-
    rbcL_species %>%
    filter(in_plastomes == TRUE) %>%
    pull(tip_name)
  
  # Exclude said species
  rbcL_sub <- pterido_rbcl_clean_seqs[!names(pterido_rbcl_clean_seqs) %in% rbcL_to_exclude]
  
  # Rename rbcL sequences as species
  rbcL_species_sub <- 
    tibble(tip_name = names(rbcL_sub)) %>%
    left_join(rbcL_species, by = "tip_name")
  
  names(rbcL_sub) <- rbcL_species_sub$species
  
  # Combine plastome and GenBank individual rbcL sequences
  rbcL_combined <- c(plastome_rbcL, rbcL_sub)
  
  # Make final tibble with all accession numbers
  rbcL_combined_accessions <- bind_rows(
    transmute(plastome_rbcL_names, species, accession, from_plastome = TRUE),
    transmute(rbcL_species_sub, species, accession, from_plastome = FALSE)
  ) %>%
    # Double check that all species are present and unique
    assert(not_na, species) %>%
    assert(is_uniq, species)
  
  # Make sure that all species are accounted for in both results
  assertthat::assert_that(
    isTRUE(
      all.equal(
        rbcL_combined_accessions %>% pull(species) %>% sort,
        rbcL_combined %>% names %>% sort
      )))
  
  return(
    list(
      rbcL_combined = rbcL_combined,
      rbcL_combined_accessions = rbcL_combined_accessions
    )
  )
  
}

#' Rename a single alignment from accessions to species names
#'
#' @param alignment Alignment to rename; matrix of class "DNAbin"
#' @param name_metadata Tibble with columns "accession" matching the
#' accessions in the alignment and "species", which will be used for renaming
#'
#' @return Matrix of class "DNAbin"
#' 
rename_alignment <- function(alignment, name_metadata) {
  
  # Only use the unique set of accessions and species names
  # (no repeats, or will mess up join.)
  name_metadata <-
    name_metadata %>%
    select(accession, species) %>%
    unique %>%
    assert(not_na, accession, species) %>%
    assert(is_uniq, accession)
  
  lookup_tibble <-
    tibble(accession = rownames(alignment)) %>%
    # Remove '_R_' that may be added to sequence names by MAFFT
    mutate(accession = str_remove_all(accession, "_R_")) %>%
    # All accessions within each gene should be unique
    assert(is_uniq, accession) %>%
    left_join(
      select(name_metadata, accession, species),
      by = "accession") %>%
    # Replace spaces with underscores
    mutate(species = str_replace_all(species, " ", "_")) %>%
    assert(not_na, species)
  
  rownames(alignment) <- lookup_tibble$species
  
  alignment
}


#' Rename a list of alignments from accessions to species names
#'
#' @param alignment_list List of alignment to rename; each must be matrix of class "DNAbin"
#' @param name_metadata Tibble with columns "accession" matching the
#' accessions in the alignment and "species", which will be used for renaming
#'
#' @return List
rename_alignment_list <- function(alignment_list, name_metadata) {
  
  purrr::map(
    alignment_list, rename_alignment,
    name_metadata = name_metadata)
  
}

#' Remove the "_R_" from mafft-generated alignments
#'
#' @param matrix Alignment matrix
#'
#' @return matrix
#' 
remove_mafft_r <- function (matrix) {
  rownames(matrix) <- rownames(matrix) %>% str_remove_all("_R_")
  matrix
}

#' Concatenate rbcL with other plastid genes
#'
#' @param plastid_genes_aligned_trimmed_renamed List of aligned genes from plastomes
#' @param rbcL_combined_aligned rbcL alignment including both GenBank individual
#' rbcL sequences and rbcL sequences from plastomes
#'
#' @return Matrix of class "DNAbin"
#' 
concatenate_rbcL_with_other_plastid_genes <- function (plastid_genes_aligned_trimmed_renamed, rbcL_combined_aligned) {
  
  # Remove plastome-only rbcL from plastome genes list
  plastid_genes_no_rbcL <- plastid_genes_aligned_trimmed_renamed %>%
    magrittr::extract(names(.) != "rbcL")
  
  # Combine combined genbank and plastome rbcL with other plastome genes
  final_gene_list <- c(
    list(rbcL = rbcL_combined_aligned),
    plastid_genes_no_rbcL)
  
  concatenate_genes(final_gene_list)
  
}

# Dating with treePL ----

#' Read in calibration and configure dates for treepl
#'
#' @param date_file_path Path to CSV file with treepl dates.
#' Must include at least the following columns:
#' - clade: name of clade
#' - taxon_1: representative taxon 1
#' - taxon_2: representative taxon 2. The MRCA of the two taxa defines the clade
#' - age: Age to assign to clade (in millions of years)
#' - age_type: 'min', 'max' or 'fixed'
#' 
#' @return Tibble with columns for use in treepl config file.
#'
load_calibration_dates <- function(date_file_path) {
  read_csv(date_file_path) %>%
    janitor::clean_names() %>%
    select(clade, age, age_type, taxon_1, taxon_2) %>%
    assert(is_uniq, clade) %>%
    assert(not_na, clade, age, age_type, taxon_1, taxon_2) %>%
    mutate(mrca = glue("mrca = {clade} {taxon_1} {taxon_2}")) %>%
    mutate(
      min_dates = case_when(
        age_type %in% c("min", "fixed") ~ glue("min = {clade} {age}"),
        TRUE ~ NA_character_),
      max_dates = case_when(
        age_type %in% c("max", "fixed") ~ glue("max = {clade} {age}"),
        TRUE ~ NA_character_))
} 

#' Do an initial treepl run to determine optimal
#' smoothing parameters with random cross-validation.
#' 
#' Requires treepl to be installed and on PATH.
#' 
#' For more details on treepl options, see
#' https://github.com/blackrim/treePL/wiki
#'
#' @param phy List of class "phy"; phylogeny
#' @param alignment List of class "DNAbin"; alignment
#' @param calibration_dates Calibration points read in with
#' load_calibration_dates()
#' @param write_tree Logical; should the phylogeny be written to working
#' directory before running treePL?
#' @param cvstart Start value for cross-validation
#' @param cvstop Stop value for cross-validation
#' @param cvsimaniter Start the number of cross validation simulated annealing 
#' iterations, default = 5000 for cross-validation
#' @param plsimaniter the number of penalized likelihood simulated annealing 
#' iterations, default = 5000
#' @param seed Seed for random number generator
#' @param thorough Logical; should the "thorough" setting in
#' treePL be used?
#' @param wd Working directory to run all treepl analyses
#' @param echo Logical; should the output be printed to the screen?
#'
run_treepl_cv <- function (
  phy, alignment, calibration_dates, 
  write_tree = TRUE,
  cvstart = "1000", cvstop = "0.1",
  cvsimaniter = "5000", 
  plsimaniter = "5000",
  nthreads = "1",
  seed,
  thorough = TRUE, wd, echo) {
  
  # Check that all taxa are in tree
  taxa <- c(calibration_dates$taxon_1, calibration_dates$taxon_2) %>% unique
  
  assertthat::assert_that(all(taxa %in% phy$tip.label),
                          msg = glue(
                            "Taxa in calibration dates not present in tree: 
                            {taxa[!taxa %in% phy$tip.label]}"))
  
  # Write tree to working directory
  phy_name <- deparse(substitute(phy))
  phy_path <- fs::path_ext_set(phy_name, "tre")
  if(isTRUE(write_tree)) {ape::write.tree(phy, fs::path(wd, phy_path))}
  
  # Get number of sites in alignment
  num_sites <- dim(alignment)[2]
  
  outfile_path <- fs::path_ext_set(paste0(phy_name, "_cv"), "out")
  
  # Write config file to working directory
  treepl_config <- c(
    glue("treefile = {phy_path}"),
    glue("numsites = {num_sites}"),
    calibration_dates$mrca,
    calibration_dates %>% filter(!is.na(min_dates)) %>% pull(min_dates),
    calibration_dates %>% filter(!is.na(max_dates)) %>% pull(max_dates),
    glue("cvstart = {cvstart}"),
    glue("cvstop = {cvstop}"),
    glue("cvsimaniter = {cvsimaniter}"),
    glue("plsimaniter = {plsimaniter}"),
    glue("cvoutfile = {outfile_path}"),
    glue("seed = {seed}"),
    glue("nthreads = {nthreads}"),
    "randomcv"
  )
  
  if(thorough) treepl_config <- c(treepl_config, "thorough")
  
  config_file_name <- glue::glue("{phy_name}_treepl_cv_config")
  
  readr::write_lines(treepl_config, fs::path(wd, config_file_name))
  
  # Run treePL
  processx::run("treePL", config_file_name, wd = wd, echo = echo)
  
  # Return cross-validation results
  read_lines(fs::path(wd, outfile_path))
  
}

#' Do an initial treepl run to determine optimal
#' smoothing parameters with random cross-validation.
#' 
#' Requires treepl to be installed and on PATH.
#' 
#' For more details on treepl options, see
#' https://github.com/blackrim/treePL/wiki
#'
#' @param phy List of class "phy"; phylogeny
#' @param alignment List of class "DNAbin"; alignment
#' @param calibration_dates Calibration points read in with
#' load_calibration_dates()
#' @param write_tree Logical; should the phylogeny be written to working
#' directory before running treePL? Can be FALSE if it is already there from
#' run_treepl_initial().
#' @param cv_results Output of run_treepl_cv() so the best smoothing
#' parameter can be selected.
#' @param cvsimaniter Start the number of cross validation simulated annealing 
#' iterations, default = 5000 for cross-validation
#' @param plsimaniter the number of penalized likelihood simulated annealing 
#' iterations, default = 5000
#' @param seed Seed for random number generator
#' @param thorough Logical; should the "thorough" setting in
#' treePL be used?
#' @param wd Working directory to run all treepl analyses
#' @param echo Logical; should the output be printed to the screen?
#'
run_treepl_prime <- function (
  phy, alignment, calibration_dates, 
  cv_results,
  write_tree = FALSE,
  cvsimaniter = "5000", 
  plsimaniter = "5000",
  nthreads = "1",
  seed,
  thorough = TRUE, wd, echo) {
  
  # Check that all taxa are in tree
  taxa <- c(calibration_dates$taxon_1, calibration_dates$taxon_2) %>% unique
  
  assertthat::assert_that(all(taxa %in% phy$tip.label),
                          msg = glue(
                            "Taxa in calibration dates not present in tree: 
                            {taxa[!taxa %in% phy$tip.label]}"))
  
  # Write tree to working directory
  phy_name <- deparse(substitute(phy))
  phy_path <- fs::path_ext_set(phy_name, "tre")
  if(isTRUE(write_tree)) {ape::write.tree(phy, fs::path(wd, phy_path))}
  
  # Get number of sites in alignment
  num_sites <- dim(alignment)[2]
  
  # Get best smoothing parameter
  # Select the optimum smoothing value (smallest chisq) from cross-validation
  # The raw output looks like this:
  # chisq: (1000) 6.7037e+30
  # chisq: (100) 3673.45
  # etc.
  best_smooth <-
    tibble(cv_result = cv_results) %>%
    mutate(smooth = str_match(cv_result, "\\((.*)\\)") %>% 
             magrittr::extract(,2) %>%
             parse_number()) %>%
    mutate(chisq = str_match(cv_result, "\\) (.*)$") %>% 
             magrittr::extract(,2) %>%
             parse_number()) %>%
    arrange(chisq) %>%
    slice(1) %>%
    pull(smooth)
  
  # Write config file to working directory
  treepl_config <- c(
    glue("treefile = {phy_path}"),
    glue("numsites = {num_sites}"),
    calibration_dates$mrca,
    calibration_dates %>% filter(!is.na(min_dates)) %>% pull(min_dates),
    calibration_dates %>% filter(!is.na(max_dates)) %>% pull(max_dates),
    glue("cvsimaniter = {cvsimaniter}"),
    glue("plsimaniter = {plsimaniter}"),
    glue("seed = {seed}"),
    glue("smooth = {best_smooth}"),
    glue("nthreads = {nthreads}"),
    "prime"
  )
  
  if(thorough) treepl_config <- c(treepl_config, "thorough")
  
  config_file_name <- glue::glue("{phy_name}_treepl_prime_config")
  
  readr::write_lines(treepl_config, fs::path(wd, config_file_name))
  
  # Run treePL
  results <- processx::run("treePL", config_file_name, wd = wd, echo = echo)
  
  # Return stdout
  read_lines(results$stdout)
  
}

#' Run treePL
#' 
#' Requires treepl to be installed and on PATH.
#' 
#' For more details on treepl options, see
#' https://github.com/blackrim/treePL/wiki
#'
#' @param phy List of class "phy"; phylogeny
#' @param alignment List of class "DNAbin"; alignment
#' @param calibration_dates Calibration points read in with
#' load_calibration_dates()
#' @param priming_results Results of running treePL with `prime` option to
#' determine optional parameters. Output of run_treepl_prime().
#' @param cv_results Results of running treePL with cross-validation to
#' determine optimal rate-smoothing parameter. Output of run_treepl_cv().
#' @param write_tree Logical; should the phylogeny be written to working
#' directory before running treePL? Can be FALSE if it is already there from
#' run_treepl_initial().
#' @param cvsimaniter Start the number of cross validation simulated annealing 
#' iterations, default = 5000 for cross-validation
#' @param plsimaniter the number of penalized likelihood simulated annealing 
#' iterations, default = 5000
#' @param seed Seed for random number generator
#' @param thorough Logical; should the "thorough" setting in
#' treePL be used?
#' @param wd Working directory to run all treepl analyses
#' @param echo Logical; should the output be printed to the screen?
#'
run_treepl <- function (
  phy, alignment, calibration_dates, 
  priming_results,
  cv_results,
  write_tree = FALSE,
  cvsimaniter = "5000", 
  plsimaniter = "5000",
  nthreads = "1",
  seed,
  thorough = TRUE, wd, echo) {
  
  # Check that all taxa are in tree
  taxa <- c(calibration_dates$taxon_1, calibration_dates$taxon_2) %>% unique
  
  assertthat::assert_that(all(taxa %in% phy$tip.label),
                          msg = glue(
                            "Taxa in calibration dates not present in tree: 
                            {taxa[!taxa %in% phy$tip.label]}"))
  
  # Write tree to working directory
  phy_name <- deparse(substitute(phy))
  phy_path <- fs::path_ext_set(phy_name, "tre")
  if(isTRUE(write_tree)) {ape::write.tree(phy, fs::path(wd, phy_path))}
  
  # Get number of sites in alignment
  num_sites <- dim(alignment)[2]
  
  # Get best smoothing parameter
  # Select the optimum smoothing value (smallest chisq) from cross-validation
  # The raw output looks like this:
  # chisq: (1000) 6.7037e+30
  # chisq: (100) 3673.45
  # etc.
  best_smooth <-
    tibble(cv_result = cv_results) %>%
    mutate(smooth = str_match(cv_result, "\\((.*)\\)") %>% 
             magrittr::extract(,2) %>%
             parse_number()) %>%
    mutate(chisq = str_match(cv_result, "\\) (.*)$") %>% 
             magrittr::extract(,2) %>%
             parse_number()) %>%
    arrange(chisq) %>%
    slice(1) %>%
    pull(smooth)
  
  # Set name of output file
  outfile_path <- fs::path_ext_set(paste0(phy_name, "_dated"), "tre")
  
  # Write config file to working directory
  treepl_config <- c(
    glue("treefile = {phy_path}"),
    glue("numsites = {num_sites}"),
    glue("smooth = {best_smooth}"),
    calibration_dates$mrca,
    calibration_dates %>% filter(!is.na(min_dates)) %>% pull(min_dates),
    calibration_dates %>% filter(!is.na(max_dates)) %>% pull(max_dates),
    glue("cvsimaniter = {cvsimaniter}"),
    glue("plsimaniter = {plsimaniter}"),
    glue("seed = {seed}"),
    glue("nthreads = {nthreads}"),
    glue("outfile = {outfile_path}"),
    # Include specifications from priming
    priming_results %>% magrittr::extract(., str_detect(., "opt =")),
    priming_results %>% magrittr::extract(., str_detect(., "optad =")),
    priming_results %>% magrittr::extract(., str_detect(., "optcvad ="))
  )
  
  if(thorough) treepl_config <- c(treepl_config, "thorough")
  
  config_file_name <- glue::glue("{phy_name}_treepl_config")
  
  readr::write_lines(treepl_config, fs::path(wd, config_file_name))
  
  # Run treePL
  processx::run("treePL", config_file_name, wd = wd, echo = echo)
  
  # Read in tree
  ape::read.tree(fs::path(wd, outfile_path))
  
}

#' Remove duplicate sequences from an alignment or tree,
#' and replace with representatives of sequence groups
#' (groups of identical sequences)
#'
#' @param rbcL_combined_alignment Combined rbcL alignment including
#' all rbcL and platome sequences
#' @param plastome_metadata_renamed Source of names of plastome sequences 
#' @param calibration_dates Tibble of taxa to be used for dating.
#' Need to check that calibration taxa are not amongst those that
#' will get removed by de-duplication.
#'
#' @return List, including:
#' grouped_alignment = Alignment with duplicate seqs removed and relabeled as groups
#' group_table = Tibble matching representative sequences to original duplicate
#' sequences so they can be added back in later.
remove_dup_seqs <- function (
  rbcL_combined_alignment, 
  plastome_metadata_renamed,
  calibration_dates) {
  
  ### Cluster identical sequences ###
  
  # It takes a very long time to run distance matrix on the full
  # alignment (> 50,000 characters), and I'm pretty sure we don't
  # have any identical plastome sequences. So remove these from
  # rbcL_combined_alignment, and check for identical seqs in the
  # remainder.
  
  # Replace spaces with underscores in species names
  # for matching to alignment names
  plastome_species <- plastome_metadata_renamed$species %>% str_replace_all(" ", "_") %>% unique
  
  rbcL_to_check <- rbcL_combined_alignment[
    !rownames(rbcL_combined_alignment) %in% plastome_species,]
  
  # Calculate raw distance matrix with pairwise deletion
  # of sites with missing data
  dist_mat <- dist.dna(rbcL_to_check, pairwise.deletion = TRUE, model = "raw")
  
  # Cluster the distances
  clusters <- hclust(dist_mat)
  
  # Find the minimum non-zero height
  min_diff_height <-
    clusters$height %>% unique %>% sort %>% magrittr::extract(2)
  
  # Find groups of sequences with minimum non-zero height.
  # i.e., groups of identical sequences
  seq_groups <- cutree(clusters, h = min_diff_height)
  
  # Make a tibble mapping each taxon to its sequence group,
  # only keep those with > 1 taxon (duplicate seqs)
  seq_groups_tibble <- 
    tibble(seq_group = seq_groups) %>%
    mutate(taxon = names(seq_groups)) %>%
    mutate(label = glue("group_{seq_group}")) %>%
    add_count(seq_group) %>%
    filter(n > 1) %>%
    # Make sure none of the duplicate sequences are used for date calibration
    assert(function(x) !x %in% calibration_dates$taxon_1, taxon) %>%
    assert(function(x) !x %in% calibration_dates$taxon_2, taxon) %>%
    arrange(seq_group, taxon) %>%
    mutate(rep_seq = ifelse(duplicated(seq_group), FALSE, TRUE))
  
  # Make a tibble of "representative" taxa, one per seq group
  # (arbitrarily choose by alph. order)
  rep_taxa_tibble <- seq_groups_tibble %>%
    filter(rep_seq == TRUE)
  
  ### Remove dups from alignment ###
  
  # Subset alignment to only "representative" seqs, rename
  # by sequence group
  pterido_rbcl_aln_rep_taxa <-
    rbcL_combined_alignment[rep_taxa_tibble$taxon, ]
  
  rownames(pterido_rbcl_aln_rep_taxa) <- rep_taxa_tibble$label
  
  # Remove all dup. seqs from alignment
  pterido_rbcl_aln_dups_removed <-
    rbcL_combined_alignment[!rownames(rbcL_combined_alignment) %in% seq_groups_tibble$taxon, ]
  
  # Make combined alignment with uniq seqs plus "rep" seqs named by seq group
  rbcL_combined_alignment_grouped <-
    rbind(pterido_rbcl_aln_dups_removed, pterido_rbcl_aln_rep_taxa)
  
  # Double-check: make sure all sequences are present in alignment
  # after restoring dup taxa
  tibble(label = rownames(rbcL_combined_alignment_grouped)) %>%
    left_join(seq_groups_tibble) %>%
    mutate(taxon = ifelse(is.na(taxon), label, taxon)) %>%
    assert(not_na, taxon) %>%
    assert(is_uniq, taxon) %>%
    verify(isTRUE(
      all.equal(
        sort(taxon),
        sort(rownames(rbcL_combined_alignment))
      )
    ), success_fun = success_logical)
  
  ### Format output ###
  # Return list: de-duplicated alignment and 
  # tibble matching group names to original sequences to add them
  # back in later
  list(
    grouped_alignment = rbcL_combined_alignment_grouped,
    group_table = seq_groups_tibble
  )
  
}

# Plotting ----

# Hacked version of ape::add.scale.bar() that allows adding units
add_scale_bar <- function (x, y, length = NULL, ask = FALSE, lwd = 1, lcol = "black", units = NA_character_,
                           ...) 
{
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  direc <- lastPP$direction
  if (is.null(length)) {
    nb.digit <- if (direc %in% c("rightwards", "leftwards")) 
      diff(range(lastPP$xx))
    else diff(range(lastPP$yy))
    length <- pretty(c(0, nb.digit)/6, 1)[2]
  }
  if (ask) {
    cat("\nClick where you want to draw the bar\n")
    x <- unlist(locator(1))
    y <- x[2]
    x <- x[1]
  }
  else if (missing(x) || missing(y)) {
    if (lastPP$type %in% c("phylogram", "cladogram")) {
      switch(direc, rightwards = {
        x <- 0
        y <- 1
      }, leftwards = {
        x <- max(lastPP$xx)
        y <- 1
      }, upwards = {
        x <- max(lastPP$xx)
        y <- 0
      }, downwards = {
        x <- 1
        y <- max(lastPP$yy)
      })
    }
    else {
      direc <- "rightwards"
      x <- lastPP$x.lim[1]
      y <- lastPP$y.lim[1]
    }
  }
  
  label <- jntools::paste3(as.character(length), units, sep = " ")
  
  switch(direc, rightwards = {
    segments(x, y, x + length, y, col = lcol, lwd = lwd)
    text(x + length * 1.1, y, label, adj = c(0, 
                                             0.5), ...)
  }, leftwards = {
    segments(x - length, y, x, y, col = lcol, lwd = lwd)
    text(x - length * 1.1, y, label, adj = c(1, 
                                             0.5), ...)
  }, upwards = {
    segments(x, y, x, y + length, col = lcol, lwd = lwd)
    text(x, y + length * 1.1, label, adj = c(0, 
                                             0.5), srt = 90, ...)
  }, downwards = {
    segments(x, y - length, x, y, col = lcol, lwd = lwd)
    text(x, y - length * 1.1, label, adj = c(0, 
                                             0.5), srt = 270, ...)
  })
}

#' Get most recent common ancestor for a selected clade
#'
#' @param phy Phylogeny
#' @param tips Dataframe with tip names and clades
#' @param clade_select Selected clade
#'
#' @return Single number
#' 
get_clade_mrca <- function (phy, tips, clade_select) {
  ape::getMRCA(
    phy, 
    tips %>% filter(clade == clade_select) %>% pull(tip)
  )
}

# Trimmomatic ----

#' Get vector of all raw fastq files for trimmomatic
#' 
#' Paths are relative to here()
#'
#' @param fastq_raw_dir Path to directory containing raw fastq files
#' from sequencing run.
#' @param fastq_pattern Grep pattern used to identify fastq files within
#' `fastq_raw_dir`
#'
#' @return vector
#' @example
#' get_fastq_files_for_trimmomatic (
#' fastq_raw_dir = "data_raw/miseq_v3_2018-12-19/indexed_reads/",
#' fastq_pattern = "fastq.gz"
#' )
get_fastq_files_for_trimmomatic <- function (fastq_raw_dir, fastq_pattern = "fastq") {
  
  list.files(fastq_raw_dir,
             full.names = TRUE,
             recursive = TRUE,
             pattern = fastq_pattern) %>%
    fs::path_abs()
  
}

#' Make a tibble of arguments to pass to trimmomatic wrapper
#' 
#' Only works on paired fastq reads differentiated by 'R1' and 'R2'.
#'
#' @param fastq_files Vector of fastq files to use as input. Paths must
#' be relative to working directory, AKA `here()`.
#' @param ids Vector of genomic IDs matching part of fastq file names
#' @param outdir Folder to write output files to.
#'
#' @return tibble
#' 
make_trimmomatic_arg_table <- function (fastq_files, ids, outdir = NULL) {
  
  assertthat::assert_that(is.character(fastq_files))
  assertthat::assert_that(is.character(ids))
  assertthat::assert_that(is.character(outdir) | is.null(outdir))
  
  extract_name <- function (name, ids) {
    name %>% str_match(paste(ids, collapse = "|")) %>% magrittr::extract(1)
  }
  
  table <-
    tibble(
      input_forward = fastq_files %>% magrittr::extract(str_detect(., "_R1_")) %>% sort,
      input_reverse = fastq_files %>% magrittr::extract(str_detect(., "_R2_")) %>% sort
    ) %>%
    mutate(
      genomicID_forward = input_forward %>% fs::path_file() %>% map_chr(~extract_name(., ids = ids)),
      genomicID_reverse = input_reverse %>% fs::path_file() %>% map_chr(~extract_name(., ids = ids))
    ) %>%
    verify(genomicID_forward == genomicID_reverse) %>%
    assert(is_uniq, genomicID_forward, genomicID_reverse) %>%
    assert(not_na, genomicID_forward, genomicID_reverse) %>%
    select(genomicID = genomicID_forward, starts_with("input")) %>%
    assert(is_uniq, genomicID) %>%
    mutate(
      output_paired_forward = glue::glue("{genomicID}_R1.fastq"),
      output_unpaired_forward = glue::glue("{genomicID}_R1_unpaired.fastq"),
      output_paired_reverse = glue::glue("{genomicID}_R2.fastq"),
      output_unpaired_reverse = glue::glue("{genomicID}_R2_unpaired.fastq")
    ) %>%
    select(
      inputFile1 = input_forward,
      inputFile2 = input_reverse,
      outputFile1P = output_paired_forward,
      outputFile1U = output_unpaired_forward,
      outputFile2P = output_paired_reverse,
      outputFile2U = output_unpaired_reverse
    )
  
  if(!is.null(outdir)) {
    table <- 
      table %>%
      mutate_at(vars(outputFile1P, outputFile1U, outputFile2P, outputFile2U), ~fs::path(outdir, .))
  }
  
  return(table)
  
}

#' Map trimmomatic over a list of fastq files.
#' 
#' Only works on paired-end fastq reads.
#'
#' @param trimmomatic_argument_table A tibble with the following columns corresponding
#' to the input to `trimmomatic_pe`: "inputFile1", "inputFile2", "outputFile1P",
#' "outputFile1U", "outputFile2P", "outputFile2U". Output of `make_trimmomatic_arg_table`.
#' @param trim_settings String with trim settigs formatted for trimmomatic, e.g.,
#' "ILLUMINACLIP:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50".
#' @param parallel Logical; should the loop be run in parallel with `furrr`? Note
#' that type of parallelization may need to be set with `furrr::plan()` first.
#' @param parallel_threads Number of threads to use per process if running in parallel.
#' @param single_threads Number of threads to use if running on one set of reads at
#' a time.
#' @param wd Working directory containing fastq files (may be in subdirectories within wd).
#' Must be an absolute path, starting with a slash.
#' @param echo Logical; should STDOUT and STERR be printed?
#'
#' @return A tibble of stdout and stderr from each call to trimmomatic.
map_trimmomatic <- function (trimmomatic_argument_table, 
                             trim_settings,
                             parallel,
                             parallel_threads,
                             single_threads,
                             wd,
                             echo) {
  if (isTRUE(parallel)) {
    
    furrr::future_pmap_dfr(
      trimmomatic_argument_table, 
      trimmomatic_pe,
      wd = wd,
      trim_settings = trim_settings,
      threads = parallel_threads,
      echo = echo)
    
  } else {
    
    pmap_dfr(
      trimmomatic_argument_table, 
      trimmomatic_pe,
      wd = wd,
      trim_settings = trim_settings,
      threads = single_threads,
      echo = echo)
    
  }
}

# Helper function to extract useful metadata from trimmomatic stderr
extract_trimmo_metadata <- function(stderr, metadata_name) {
  str_extract(stderr, glue("{metadata_name}: [[:digit:]]*")) %>% 
    str_remove(glue("{metadata_name}: ")) %>% 
    as.numeric
}

# Tidy raw trimmomatic results from map_trimmomatic
tidy_trimmomatic <- function (trimmomatic_results_raw) { 
  trimmomatic_results_raw %>%
    mutate(
      input_read_pairs = extract_trimmo_metadata(stderr, "Input Read Pairs"),
      both_surviving = extract_trimmo_metadata(stderr, "Both Surviving"),
      forward_only_surviving = extract_trimmo_metadata(stderr, "Forward Only Surviving"),
      reverse_only_surviving = extract_trimmo_metadata(stderr, "Reverse Only Surviving"),
      dropped = extract_trimmo_metadata(stderr, "Dropped"),
      sample = str_extract(stderr, "\\S*_R1.fastq") %>% str_remove("_R1.fastq") %>% fs::path_file() %>% as.character
    ) %>% 
    select(-stderr) %>%
    mutate(both_surviving_frac = both_surviving / input_read_pairs) %>%
    select(sample, everything())
}

#' Run trimmomatic in paired-end mode.
#' 
#' Requires either trimmomatic to be installed and in the user's PATH, or 
#' docker to be installed.
#' 
#' @param threads Numeric; specify the number of threads to use (optional).
#' @param phred_type String; specify the base quality encoding. Either "phred33" (default) or "phred64".
#' @param trimLogFile String; name to write log file (optional).
#' @param statsSummaryFile String; name to write summary stats file (optional).
#' @param quiet Logical; should trimmomatic be run in quiet mode?
#' @param validatePairs Logical: should trimmomatic validate pairs?
#' @param inputFile1 String; file name of input forward reads.
#' @param inputFile2 String; file name of input reverse reads.
#' @param outputFile1P String; file name of output trimmed paired forward reads.
#' @param outputFile1U String; file name of output trimmed unpaired forward reads.
#' @param outputFile2P String; file name of output trimmed paired reverse reads.
#' @param outputFile2U String; file name of output trimmed unpaired forward reads.
#' @param trim_settings String; trim settings (e.g, 
#' "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36").
#' For details on specifying trim settings, see trimmomatic manual. 
#' @param ... Other arguments; not used by this function, but meant for tracking by drake.
#'
#' @return A list with components specified in {\link[processx]{run}. Externally, trimmed reads
#' will be written in the working directory.
#' @examples
#' \dontrun{
#' trimmomatic_pe(wd = here::here(),
#' inputFile1 = "input_forward.fq.gz",
#' inputFile2 = "input_reverse.fq.gz",
#' outputFile1P = "output_forward_paired.fq.gz",
#' outputFile1U = "output_forward_unpaired.fq.gz",
#' outputFile2P = "output_reverse_paired.fq.gz",
#' outputFile2U = "output_reverse_unpaired.fq.gz",
#' trim_settings = "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
#' }
trimmomatic_pe <- function (threads = NULL, 
                            phred_type = "phred33", 
                            trimLogFile = NULL, 
                            statsSummaryFile = NULL, 
                            quiet = FALSE,
                            validatePairs = FALSE,
                            inputFile1,
                            inputFile2,
                            outputFile1P,
                            outputFile1U,
                            outputFile2P,
                            outputFile2U,
                            trim_settings, ...) {
  
  # Make sure input types are correct
  if(!is.null(threads)) assertthat::assert_that(assertthat::is.number(threads))
  if(!is.null(trimLogFile)) assertthat::assert_that(assertthat::is.string(trimLogFile))
  if(!is.null(statsSummaryFile)) assertthat::assert_that(assertthat::is.string(statsSummaryFile))
  assertthat::assert_that(is.logical(quiet))
  assertthat::assert_that(is.logical(validatePairs))
  assertthat::assert_that(assertthat::is.string(inputFile1))
  assertthat::assert_that(assertthat::is.string(inputFile2))
  assertthat::assert_that(assertthat::is.string(outputFile1P))
  assertthat::assert_that(assertthat::is.string(outputFile1U))
  assertthat::assert_that(assertthat::is.string(outputFile2P))
  assertthat::assert_that(assertthat::is.string(outputFile2U))
  assertthat::assert_that(assertthat::is.string(trim_settings))
  
  # Modify arguments for processx::run()
  threads <- if(is.null(threads)) NULL else glue::glue("-threads {threads}")
  phred_type <- if(is.null(phred_type)) NULL else glue::glue("-{phred_type}")
  trimLogFile <- if(is.null(trimLogFile)) NULL else glue::glue("-trimlog {trimLogFile}")
  statsSummaryFile <- if(is.null(statsSummaryFile)) NULL else glue::glue("-summary {statsSummaryFile}")
  quiet <- if(!isTRUE(quiet)) NULL else "-quiet"
  validatePairs <- if(!isTRUE(validatePairs)) NULL else "-validatePairs"
  
  trimmomatic_arguments <- c("PE", threads, phred_type, trimLogFile, 
                             inputFile1, inputFile2, 
                             outputFile1P, outputFile1U,
                             outputFile2P, outputFile2U,
                             trim_settings)
  
  # Return stdout to character vector, then parse this into a tibble
  output <- system2("trimmomatic", trimmomatic_arguments, stdout = TRUE, stderr = TRUE)
  
  tibble::tibble(stderr = unlist(output) %>% paste(collapse = " ")) %>%
    tidy_trimmomatic
}

# Hypbiper ----

#' Get taxon name for a genbank accession number
#'
#' @param accession String; genbank accession number
#' @param rank String; taxonomic level. Must be one of
#' "superkingdom", "kingdom", "phylum", "subphylum",
#' "class", "subclass", "order", "suborder"
#' "family", "genus", "species"
#'
#' @return String
#' @examples
#' get_taxon_from_gb("KY427331") 
get_taxon_from_gb <- function (accession, rank = "species") {
  
  # Sleep for 1/3 second to avoid triggering Error: Too Many Requests (HTTP 429)
  # May have something to do with this?
  # https://github.com/ropensci/taxize/issues/722
  Sys.sleep(0.33)
  
  assertthat::assert_that(assertthat::is.string(accession))
  assertthat::assert_that(assertthat::is.string(rank))
  
  # Check for ENTREZ KEY
  assertthat::validate_that(
    nchar(Sys.getenv("ENTREZ_KEY")) > 0,
    msg = "Warning: ENTREZ_KEY not detected. See taxize::key_helpers()."
  )
  
  taxonomy <-
    taxize::genbank2uid(id = accession) %>%
    taxize::classification(db = "ncbi") %>%
    purrr::flatten()
  
  assertthat::assert_that(rank %in% taxonomy$rank)
  taxonomy$name[taxonomy$rank == rank]
}

# Count the number of missing bases in a DNA sequence
# Here "missing" is anything other than A,C,T, or G
#' @param seq List of class DNAbin of length one.
count_missing <- function (seq) {
  seq <- as.character(seq)
  non_missing <- sum(str_count(seq, "a|t|c|g|A|T|C|G"))
  total_chars <- sum(nchar(seq))
  total_chars - non_missing
}

#' Construct plastid HybPiper target file
#'
#' @param gene_list Named list of genes, each which contains a list of DNA sequences
#' of class DNAbin. Output of assemble_gene_set().
#' @param accessions Vector of accessions.
#' @param gene_names Vector of gene names.
#' @param taxonomy_data Dataframe of taxonomic data for ferns following PPGI
#' system. Hosted at http://bit.ly/ppgi_taxonomy
#' @param out_path Path to write out HybPiper target file.
#'
#' @return List of class DNAbin
#'
#' @examples
#' test_accessions <- c("KP136830", "KY427346")
#' test_genes <- c("accD", "atpA", "psbA")
#' test_gene_list <- assemble_gene_set(test_accessions, test_genes)
#' ppgi_taxonomy <- read_csv("http://bit.ly/ppgi_taxonomy")
#' make_plastid_target_file (test_gene_list, test_accessions, test_genes, ppgi_taxonomy, "test.fasta")
#' 
make_plastid_target_file <- function (gene_list, accessions, gene_names, taxonomy_data, out_path) {
  
  # Get species names for accessions, join with PPGI taxonomy.
  species <-
    accessions %>%
    set_names(.) %>%
    map_chr(get_taxon_from_gb)
  
  species_taxonomy <-
    tibble(species = species, accession = names(species)) %>%
    mutate(genus = str_split(species, " ") %>% map_chr(1)) %>%
    left_join(select(taxonomy_data, suborder, family, subfamily, genus))
  
  # Filter the list of plastid genes by taxonomy (must be a eupolypod II fern),
  # then choose one sequence per genus that is the longest seq with least missing
  # data.
  plastid_phy_gene_list <- flatten(gene_list)
  
  hybpiper_target_data <-
    tibble(seq = plastid_phy_gene_list, seq_name = names(plastid_phy_gene_list)) %>%
    separate(seq_name, into = c("accession", "gene"), sep = "-", remove = FALSE) %>%
    left_join(species_taxonomy) %>%
    filter(suborder == "Aspleniineae") %>%
    mutate(n_missing = map_dbl(seq, count_missing),
           length = map_dbl(seq, length)) %>%
    group_by(gene, genus) %>%
    filter(length == max(length)) %>%
    filter(n_missing == min(n_missing)) %>%
    slice(1) %>%
    ungroup()
  
  # These are the DNA sequences we will use as targets for HybPiper.
  plastid_targets <- pull(hybpiper_target_data, seq) %>%
    set_names(hybpiper_target_data$seq_name) %>%
    map(as.character) %>%
    ape::as.DNAbin()
  
  # Check that all names meet hybpiper formatting requirements
  # (sample and gene separated by dash).
  assertthat::assert_that(
    str_split(names(plastid_targets), "-") %>% 
      map_chr(1) %>% 
      setdiff(accessions) %>%
      length %>% 
      magrittr::equals(0),
    msg = "Target names don't match HybPiper format"
  )
  
  assertthat::assert_that(
    str_split(names(plastid_targets), "-") %>%
      map_chr(2) %>%
      setdiff(gene_names) %>%
      length %>%
      magrittr::equals(0),
    msg = "Target names don't match HybPiper format"
  )
  
  assertthat::assert_that(
    str_split(names(plastid_targets), "-") %>%
      map_dbl(length) %>%
      unique %>%
      magrittr::equals(2),
    msg = "Target names don't match HybPiper format"
  )
  
  # Write out target file for hybpiper
  ape::write.FASTA(plastid_targets, out_path)
  
  plastid_targets
  
}

# Get vector of trimmed reads for hybpiper
get_reads <- function (data_dir, pattern, ...) {
  list.files(data_dir, pattern = pattern, full.names = TRUE) %>%
    fs::path_norm() %>%
    sort
}

# Make list of readfiles for hybpiper readsfirst.py
make_paired_reads_list <- function (
  forward_reads, reverse_reads, 
  forward_read_ending = "_R1.fastq",
  reverse_read_ending = "_R2.fastq") {
  tibble(forward_reads = forward_reads) %>%
    mutate(reverse_reads = reverse_reads,
           forward_prefix = 
             fs::path_file(forward_reads) %>%
             str_remove(., forward_read_ending),
           reverse_prefix = 
             fs::path_file(reverse_reads) %>%
             str_remove(., reverse_read_ending)
    ) %>%
    # Make sure the read names are exactly the same except for the ending
    verify(forward_prefix == reverse_prefix) %>%
    mutate(readfiles = map2(forward_reads, reverse_reads, c)) %>%
    pull(readfiles)
}

#' Make list of samples analyzed by HybPiper
#' 
#' Outputs a text file consisting of sample names, one per line.
#'
#' @param in_dir Directory containing hybpiper output, including one folder per sample.
#' @param pattern Grep pattern to match hybpiper output folders.
#' @param out_path Path to write sample list.
#' @param ... Other arguments not used by this function but meant for tracking
#' with `drake`.
#' @return Character vector; externally, a text file will be written to out_path.
#' @example 
#' eupoly2_samples <- make_hybpiper_sample_file(
#'   in_dir = here("03_hybpiper"), 
#'   pattern = "4938|JNG", 
#'   out_path = file_out(here("03_hybpiper/eupoly2_samples.fasta"))
#' )
make_hybpiper_sample_file <- function (in_dir, pattern, out_path, ...) {
  list.files(in_dir, pattern = pattern) %>%
    readr::write_lines(out_path)
}

#' Run HybPiper reads_first.py
#'
#' Extracts hits to target sequences from a (pair of) input read(s). This 
#' requires HybPiper scripts to be on the user's PATH.
#'
#' @param wd String; working directory. All output will be written here.
#' @param echo Logical; should STDOUT and STERR be printed?
#' @param baitfile String; FASTA file containing bait sequences for each gene. 
#' If there are multiple baits for a gene, the id must be of the form: >Taxon-geneName.
#' (From https://github.com/mossmatters/HybPiper/blob/master/reads_first.py).
#' @param readfiles Character vector; One or more read files to start the pipeline. 
#' If exactly two are specified, will assume it is paired Illumina reads.
#' (From https://github.com/mossmatters/HybPiper/blob/master/reads_first.py).
#' @param prefix String; Directory name for pipeline output, default is to use 
#' the FASTQ file name.
#' (From https://github.com/mossmatters/HybPiper/blob/master/reads_first.py).
#' @param bwa Logical; Use BWA to search reads for hits to target. 
#' Requires BWA and a bait file that is nucleotides!
#' @param cpu Numeric; Limit the number of CPUs. Default is to use all cores available.
#' (From https://github.com/mossmatters/HybPiper/blob/master/reads_first.py).
#' @param ... Other arguments not used by this function but meant for tracking
#' with `drake`.
#'
#' @return A list with components specified in {\link[processx]{run}. 
#' Externally, results of reads_first.py will be written in the working 
#' directory.
#'
#' @examples
reads_first <- function (wd, echo = FALSE, baitfile, readfiles, prefix = NULL, bwa = FALSE, cpu = NULL, other_args = NULL, ...) {
  
  # Make sure input types are correct
  assertthat::assert_that(assertthat::is.readable(wd))
  assertthat::assert_that(is.logical(echo))
  assertthat::assert_that(assertthat::is.readable(baitfile))
  assertthat::assert_that(is.character(readfiles))
  if(!is.null(prefix)) 
    assertthat::assert_that(assertthat::is.string(prefix))
  if(!is.null(cpu)) 
    assertthat::assert_that(assertthat::is.number(cpu))
  assertthat::assert_that(is.logical(bwa))
  if(!is.null(other_args)) 
    assertthat::assert_that(is.character(other_args))
  
  # Modify arguments for processx::run()
  wd <- fs::path_abs(wd)
  baitfile <- fs::path_abs(baitfile)
  
  hybpiper_arguments <- c("-b", baitfile, 
                          "-r", readfiles,
                          if(!is.null(prefix)) "--prefix", 
                          prefix,
                          if(!is.null(cpu)) "--cpu", 
                          cpu,
                          if(isTRUE(bwa)) "--bwa",
                          other_args)
  
  # Run command
  processx::run(
    "reads_first.py", 
    hybpiper_arguments, 
    wd = wd, echo = echo)
}

#' Check the length of recovered genes by sample.
#' 
#' Runs HybPiper get_seq_lengths.py and returns results as a tibble.
#' 
#' @param echo Logical; should STDOUT and STERR be printed?
#' @param wd String; working directory. This should contain the results of HybPiper runs
#' (one folder per sample named by sample).
#' @param baitfile String; path to FASTA file containing bait sequences for each gene. 
#' If there are multiple baits for a gene, the id must be of the form: >Taxon-geneName.
#' @param namelistfile String; path to file that contains a list of all the HybPiper 
#' directories for each sample.
#' @param out_path Optional; Path to write output of get_seq_lengths.py.
#' @param sequenceType String; type of molecule used for the baitfile, either
#' "dna" or "aa".
#' @param ... Other arguments not used by this function but meant for tracking
#' with `drake`.
#'
#' @return A tibble. The first row is the mean lengths of the target genes;
#' other rows are gene lengths for each sample. If `out_path` is specified, a
#' text file will be written with the output of get_seq_lengths.py.
#'
#' @examples
#' plastid_lengths <- get_seq_lengths(
#'   baitfile = here("01_trimmomatic/old/reads_paired/eupoly2_plastid_targets.fasta"), 
#'   namelistfile = here("01_trimmomatic/old/reads_paired/samples.txt"),
#'   out_path = here("01_trimmomatic/old/reads_paired/lengths.txt"),
#'   sequenceType = "dna",
#'   wd = here("01_trimmomatic/old/reads_paired/")
#' )
#' @references https://github.com/mossmatters/HybPiper/blob/master/reads_first.py
#' 
get_seq_lengths <- function (echo = FALSE, wd, baitfile, namelistfile, sequenceType, out_path = NULL, ...) {
  
  # Make sure input types are correct
  assertthat::assert_that(is.logical(echo))
  assertthat::assert_that(assertthat::is.readable(wd))
  assertthat::assert_that(assertthat::is.readable(baitfile))
  assertthat::assert_that(assertthat::is.readable(namelistfile))
  assertthat::assert_that(
    sequenceType %in% c("dna", "aa"), 
    msg = "Must choose either 'dna' or 'aa' for sequenceType."
  )
  if(!is_null(out_path)) assertthat::assert_that(assertthat::is.string(out_path))
  
  # Modify arguments for processx::run()
  wd <- fs::path_abs(wd)
  hybpiper_arguments <- c(baitfile, 
                          namelistfile,
                          sequenceType)
  
  
  # Run command
  results <- processx::run(
    "get_seq_lengths.py", 
    hybpiper_arguments, wd = wd, echo = echo)
  
  # Check for hybpiper error
  assertthat::assert_that(
    !str_detect(results$stdout[1], "^\n\nUsage:"),
    msg = "Error by HybPiper")
  
  # Write out and convert raw text to tibble
  lengths_raw <- results$stdout
  out_path <- ifelse(is_null(out_path), tempfile(), out_path)
  # Don't use write_lines, it will add a trailing break
  readr::write_file(lengths_raw, out_path)
  readr::read_delim(out_path, delim = "\t")
  
}

#' Check HybPiper summary statistics.
#' 
#' Runs HybPiper hybpiper_stats.py and returns results as a tibble.
#' 
#' @param echo Logical; should STDOUT and STERR be printed?
#' @param wd String; working directory. This should contain the results of HybPiper runs
#' (one folder per sample named by sample).
#' @param seq_lengths String; path to text file containing lengths of each gene by sample
#' (output of get_seq_lengths).
#' @param namelistfile String; path to file that contains a list of all the HybPiper 
#' directories for each sample.
#' @param ... Other arguments not used by this function but meant for tracking
#' with `drake`.
#'
#' @return A tibble. 
#'
#' @examples
#' plastid_stats <- hybpiper_stats(
#'   seq_lengths = here("01_trimmomatic/old/reads_paired/lengths.txt"), 
#'   namelistfile = here("01_trimmomatic/old/reads_paired/samples.txt"),
#'   wd = here("01_trimmomatic/old/reads_paired/")
#' )
#' @references https://github.com/mossmatters/HybPiper/wiki/Tutorial
#' 
hybpiper_stats <- function (echo = FALSE, wd, seq_lengths, namelistfile, ...) {
  
  # Make sure input types are correct
  assertthat::assert_that(is.logical(echo))
  assertthat::assert_that(assertthat::is.readable(wd))
  assertthat::assert_that(assertthat::is.readable(seq_lengths))
  assertthat::assert_that(assertthat::is.readable(namelistfile))
  
  # Modify arguments for processx::run()
  wd <- fs::path_abs(wd)
  hybpiper_arguments <- c(seq_lengths, 
                          namelistfile)
  
  # Run command
  results <- processx::run(
    "hybpiper_stats.py", 
    hybpiper_arguments, wd = wd, echo = echo)
  
  # Check for hybpiper error
  assertthat::assert_that(
    !str_detect(results$stdout[1], "^\n\nUsage:"),
    msg = "Error by HybPiper")
  
  # Convert raw text to tibble
  lengths_raw <- results$stdout
  out_path <- tempfile()
  readr::write_file(lengths_raw, out_path)
  readr::read_delim(out_path, delim = "\t")
  
}

#' Retrieve sequences from HybPiper run
#' 
#' This script fetches the sequences recovered from the same gene 
#' for many samples and generates an unaligned multi-FASTA file for each gene.
#' 
#' @param echo Logical; should STDOUT and STERR be printed?
#' @param wd String; working directory. Fasta files will be written here.
#' @param baitfile String; path to FASTA file containing bait sequences for each gene. 
#' If there are multiple baits for a gene, the id must be of the form: >Taxon-geneName.
#' @param sequence_dir String; directory containing the results of HybPiper runs
#' (one folder per sample named by sample).
#' @param sequenceType String; type of molecule to return, either
#' "dna" or "aa". Can also choose "intron" or "supercontig" if intronerate
#' was run previously on the HybPiper runs.
#' @param ... Other arguments not used by this function but meant for tracking
#' with `drake`.
#'
#' @return A list with components specified in {\link[processx]{run}. 
#' Externally, one fasta file per gene will be written to the working directory.
#'
#' @references https://github.com/mossmatters/HybPiper/wiki/Tutorial
#' 
retrieve_sequences <- function (echo = FALSE, wd, baitfile, sequence_dir, sequenceType, ...) {
  
  # Make sure input types are correct
  assertthat::assert_that(is.logical(echo))
  assertthat::assert_that(assertthat::is.dir(wd))
  assertthat::assert_that(assertthat::is.readable(baitfile))
  assertthat::assert_that(assertthat::is.dir(sequence_dir))
  assertthat::assert_that(
    sequenceType %in% c("dna", "aa", "intron", "supercontig"), 
    msg = "Must choose from 'dna', 'aa', 'intron', or 'supercontig' for sequenceType."
  )
  
  # Modify arguments for processx::run()
  wd <- fs::path_abs(wd)
  hybpiper_arguments <- c(baitfile, 
                          sequence_dir,
                          sequenceType)
  
  # Run command
  processx::run(
    "retrieve_sequences.py", 
    hybpiper_arguments, wd = wd, echo = echo)
  
}

#' Retrieve intron sequences from HybPiper run
#' 
#' Given a completed run of reads_first.py for a sample, run this script to 
#' generate "gene" sequences for each locus. The script will generate two new 
#' sequence files for each gene:
#' - supercontig: A sequence containing all assembled contigs with a unique 
#' alignment to the reference protein, concatenated into one sequence.
#' - introns: The supercontig with the exon sequences removed.
#' 
#' @param echo Logical; should STDOUT and STERR be printed?
#' @param wd String; working directory. Fasta files will be written here.
#' @param prefix String; name of output folder from reads_first containing
#' hybpiper output.
#' @param ... Other arguments not used by this function but meant for tracking
#' with `drake`.
#'
#' @return A list with components specified in {\link[processx]{run}. 
#' Externally, one fasta file per gene will be written to the working directory.
#'
#' @references https://github.com/mossmatters/HybPiper/wiki/Tutorial
#' 
intronerate <- function (echo = FALSE, wd, prefix, ...) {
  
  # Make sure input types are correct
  assertthat::assert_that(assertthat::is.string(prefix))
  
  assertthat::assert_that(is.logical(echo))
  assertthat::assert_that(assertthat::is.dir(wd))
  
  # Modify arguments for processx::run()
  wd <- fs::path_abs(wd)
  
  hybpiper_arguments <- c("--prefix", 
                          prefix)
  
  # Run command
  processx::run(
    "intronerate.py", 
    hybpiper_arguments, wd = wd, echo = echo)
  
}

#' Get statistics on the number of reads obtained by HybPiper
#' 
#' (just reads, not assembled genes)
#'
#' @param hybpiper_dir Path to HybPiper working directory
#' @param ... Other arguments not used by this function but meant for tracking
#' with `drake`.
#'
get_read_stats <- function (hybpiper_dir, ...) {
  
  get_total_seq_len <- function (dna) {
    purrr::map_dbl(dna, length) %>% sum()
  }
  
  tibble(
    file = list.files(hybpiper_dir, pattern = "interleaved\\.fasta", recursive = TRUE, full.names = TRUE),
    data = purrr::map(file, ape::read.FASTA)
  ) %>%
    mutate(
      sample = purrr::map_chr(file, ~str_split(., "\\/") %>% map_chr(7)),
      gene = purrr::map_chr(file, ~str_split(., "\\/") %>% map_chr(8)),
      n_seq = purrr::map_dbl(data, length),
      n_bp = purrr::map_dbl(data, get_total_seq_len)
    ) %>%
    select(sample, gene, n_seq, n_bp)
  
}

