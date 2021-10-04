# Data loading ----

#' Load data from the Catalog of Life
#'
#' @param col_data_path Path to TSV downloaded from the Catalog of Life 
#' https://www.catalogueoflife.org/data/download
#'
#' @return Dataframe (tibble)
#' 
load_col <- function(col_data_path) {
  # Use data.table as it handles quotation marks in data better than read_*() functions
  data.table::fread(file = col_data_path, sep = "\t", stringsAsFactors = FALSE) %>%
    janitor::clean_names() %>%
    tibble::as_tibble()
}

#' Extract Ferns of the World taxonomic data from Catalog of Life
#' 
#' Also do some quality checks, and filter to only data at species level or below
#'
#' @param col_data Catalog of Life data. Output of of load_col()
#'
#' @return Dataframe (tibble)
#' 
extract_fow_from_col <- function(col_data) {
  # Filter Catalog of Life data to only Ferns of the World
  fow <- col_data %>% 
    filter(dwc_dataset_id == "1140") %>%
    select(
      taxon_id = dwc_taxon_id, 
      accepted_name_taxon_id = dwc_accepted_name_usage_id, 
      status = dwc_taxonomic_status,
      rank = dwc_taxon_rank, 
      scientific_name = dwc_scientific_name) %>%
    # Keep only species level and below
    filter(rank %in% c("form", "infraspecific name", "species", "subform", "subspecies", "subvariety", "variety")) %>%
    # Filter some names that were incorrectly labeled species level
    filter(str_detect(scientific_name, "Polypodiaceae tribe Thelypterideae|Asplenium grex Triblemma|Pteridaceae tribus Platyzomateae|Filicaceae tribus Taenitideae", negate = TRUE)) %>%
    mutate(status = str_replace_all(status, "provisionally accepted", "accepted")) %>%
    select(-rank) %>%
    # Note: taxon_id is unique, but scientific_name may not be (esp in case of ambiguous synonyms)
    assert(not_na, taxon_id, scientific_name) %>% 
    assert(is_uniq, taxon_id)
  
  # Make sure all synonyms map correctly
  fow_accepted <- 
    fow %>%
    filter(str_detect(status, "accepted")) %>%
    select(taxon_id, scientific_name, -status)
  
  fow_synonyms <- 
    fow %>%
    filter(str_detect(status, "synonym")) %>%
    select(taxon_id, accepted_name_taxon_id, scientific_name, -status)
  
  fow_synonyms %>%
    anti_join(fow_accepted, by = c(accepted_name_taxon_id = "taxon_id")) %>%
    verify(nrow(.) == 0, success_fun = success_logical)
  
  # Make sure all accepted names and synonyms are accounted for
  bind_rows(fow_accepted, fow_synonyms) %>%
    assert(is_uniq, taxon_id) %>%
    anti_join(fow, by = "taxon_id") %>%
    verify(nrow(.) == 0, success_fun = success_logical)
  
  # Make sure accepted names and synonyms are distinct
  # A few repeats. Leave these in for now, but will need to fix.
  # fow_accepted %>%
  #   inner_join(fow_synonyms, by = c(name = "synonym")) %>%
  #   select(scientific_name = name) %>% 
  #   left_join(fow)
  
  fow
}

# Download Sanger sequences from GenBank----

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
  
  # Check for FEATURES and ORIGIN field; if either is missing, return NULL
  if (str_detect(gb_entry, "FEATURES", negate = TRUE)) {
    message("Genbank flatfile not valid (missing FEATURES); no sequence extracted")
    return(NULL)
  }
  
  if (str_detect(gb_entry, "ORIGIN",  negate = TRUE)) {
    message("Genbank flatfile not valid (missing ORIGIN); no sequence extracted")
    return(NULL)
  }
  
  if (str_detect(gb_entry, "ACCESSION",  negate = TRUE)) {
    message("Genbank flatfile not valid (missing ACCESSION); no sequence extracted")
    return(NULL)
  }
  
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
  gene_detected <- any(str_detect(gene_range_list, regex(gene, ignore_case = TRUE)))
  gene_detected_msg <- assertthat::validate_that(
    gene_detected,
    msg = glue::glue("Gene {gene} not detected in accession {accession}")
  )
  if(!gene_detected) {
    message(gene_detected_msg)
    return(NULL)
  }
  
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
  
  # Check for duplicated genes
  gene_duplicated <- length(unique(gene_range)) <= 2
  gene_duplicated_msg <- assertthat::validate_that(
    gene_duplicated,
    msg = glue::glue("Duplicate copies of gene {gene} detected in accession {accession}")
  )
  if(!gene_duplicated) {
    message(gene_duplicated_msg)
    return(NULL)
  }
  
  # Check that full range of gene was detected
  gene_full <- length(gene_range) > 1 && 
    is.numeric(gene_range) && 
    !anyNA(gene_range) && 
    gene_range[1] <= gene_range[2] && 
    gene_range[2] >= gene_range[1]
  
  gene_full_msg <- assertthat::validate_that(
    gene_full,
    msg = glue::glue("Full range of gene {gene} not detected in accession {accession}")
  )
  
  if(!gene_full) {
    message(gene_full_msg)
    return(NULL)
  }
  
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
    # Drop any NULL values
    purrr::compact() %>%
    # Name them as the accession
    purrr::set_names(map_chr(., names)) %>%
    # Convert to ape format
    # - each sequence needs to be a character vector with each letter as an element
    map(stringr::str_split, "") %>%
    map(ape::as.DNAbin) %>%
    do.call(c, .)
  
}

#' Download a set of fern sequences for a given gene
#'
#' @param gene Name of gene
#' @param start_date Earliest date to download
#' @param end_date Most recent date to download
#' @param return_df Logical; return results as a dataframe? If FALSE, returns results
#' as a list
#'
#' @return List of class DNAbin or dataframe with one row per sequence and columns
#' for the sequence and accession
#' 
#' fetch_fern_gene("rbcL", start_date = "2018/01/01", end_date = "2018/01/10")
fetch_fern_gene <- function(gene, start_date = "1980/01/01", end_date, return_df = TRUE) {
  
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
  print(glue("Found {num_hits} sequences (UIDs) for {gene}"))
  
  # Download complete GenBank record for each and write it to a temporary file
  temp_dir <- tempdir()
  temp_file <- fs::path(temp_dir, "gb_records.txt")
  
  # Make sure temp file doesn't already exist
  if(fs::file_exists(temp_file)) fs::file_delete(temp_file)
  
  reutils::efetch(uid, "nucleotide", rettype = "gb", retmode = "text", outfile = temp_file)
  
  # Parse flatfile
  seqs <- parse_dna_from_flatfile(temp_file, gene)
  
  fs::file_delete(temp_file)
  
  # Return list 
  if(return_df == FALSE) return(seqs)
  
  # Or return dataframe
  tibble::tibble(seq = split(seqs, 1:length(seqs)), accession = names(seqs), gene = gene)
  
}

# Download metadata for Sanger sequences from GenBank ----

#' Fetch metadata from GenBank
#'
#' @param query String to use for querying GenBank
#' @param col_select Character vector; columns of metadata to retain in output
#'
#' @return Tibble
#' 
fetch_metadata <- function(
  query = NULL,
  col_select = c("gi", "caption", "taxid", "title", "slen", "subtype", "subname")) {
  
  assertthat::assert_that(assertthat::is.string(query))
  
  assertthat::assert_that(is.character(col_select))
  
  # Do an initial search without downloading any IDs to see how many hits
  # we get.
  initial_genbank_results <- rentrez::entrez_search(
    db = "nucleotide",
    term = query,
    use_history = FALSE
  )
  
  assertthat::assert_that(
    initial_genbank_results$count > 0,
    msg = "Query resulted in no hits")
  
  # Download IDs with maximum set to 1 more than the total number of hits.
  genbank_results <- rentrez::entrez_search(
    db = "nucleotide",
    term = query,
    use_history = FALSE,
    retmax = initial_genbank_results$count + 1
  )
  
  # Define internal function to download genbank data into tibble
  entrez_summary_gb <- function(id, col_select) {
    # Download data
    rentrez::entrez_summary(db = "nucleotide", id = id) %>%
      # Extract selected columns from result
      purrr::map_dfr(magrittr::extract, col_select) %>%
      # Make sure taxid column is character
      mutate(taxid = as.character(taxid)) %>%
      assert(not_na, taxid)
  }
  
  # Extract list of IDs from search results
  genbank_ids <- genbank_results$ids
  
  # Fetch metadata for each ID and extract selected columns
  if (length(genbank_ids) == 1) {
    rentrez_results <- rentrez::entrez_summary(db = "nucleotide", id = genbank_ids) %>%
      magrittr::extract(col_select) %>%
      tibble::as_tibble() %>%
      mutate(taxid = as.character(taxid)) %>%
      assert(not_na, taxid)
  } else {
    # Split input vector into chunks
    n <- length(genbank_ids)
    chunk_size <- 200
    r <- rep(1:ceiling(n/chunk_size), each=chunk_size)[1:n]
    genbank_ids_list <- split(genbank_ids, r) %>% magrittr::set_names(NULL)
    # Download results for each chunk
    rentrez_results <- map_df(genbank_ids_list, ~entrez_summary_gb(., col_select = col_select))
  }
  
  return(rentrez_results)
  
}

#' Helper function for fetch_genbank_refs() to parse character vector from GenBank record into dataframe
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
  metadata <- fetch_metadata(query) %>%
    rename(accession = caption) %>%
    # GenBank accession should be non-missing, unique
    assert(not_na, accession) %>%
    assert(is_uniq, accession)
  
  # Also fetch reference data (publication to cite for the sequence).
  # If the title is only "Direct Submission" consider this to be NA
  # (no real pub associated with the sequence, can't match on this)
  ref_data <- fetch_genbank_refs(query) %>%
    transmute(
      accession, 
      publication = str_replace_all(title, "Direct Submission", NA_character_)) %>%
    # GenBank accession should be non-missing, unique
    assert(not_na, accession) %>%
    assert(is_uniq, accession)
  
  # Combine standard metadata with publication metadata
  left_join(metadata, ref_data, by = "accession") %>%
    assert(is_uniq, accession) %>%
    mutate(gene = gene)
  
}

# Combine Sanger sequences data ----

#' Combine sanger sequence metadata with sequences, join to resolved names
#' and filter by sequence length and if name was resolved or not
#' 
#' Drops sequences with scientific names that could not be resolved, nothogenera
#'
#' @param raw_meta Sanger sequence metadata; output of fetch_fern_metadata()
#' @param raw_fasta Sanger sequences; output of fetch_fern_gene()
#' @param ncbi_accepted_names_map Dataframe mapping NCBI taxid to accepted
#' species name; output of make_ncbi_accepted_names_map()
#' @param ppgi PPGI taxonomic system
#'
#' @return Tibble with Sanger sequence metadata, sequences, and accepted name
#' 
combine_and_filter_sanger <- function(raw_meta, raw_fasta, ncbi_accepted_names_map, ppgi) {
  # Join metadata and fasta sequences
  raw_meta %>%
    left_join(raw_fasta, by = c("accession", "gene")) %>%
    # Inner join to name resolution results: will drop un-resolved names
    inner_join(ncbi_accepted_names_map, by = "taxid") %>%
    # Drop nothogenera
    mutate(genus = stringr::str_split(taxon, "_") %>% purrr::map_chr(1)) %>%
    left_join(
      select(ppgi, genus, nothogenus), 
      by = "genus") %>%
    assert(not_na, nothogenus) %>%
    filter(nothogenus == "no") %>%
    select(-genus, -nothogenus) %>%
    # Filter by minimum seq. length
    filter(slen > 400) %>%
    assert(not_na, accession, seq, accepted_name, taxon) %>%
    # Filter out null sequences
    filter(!map_lgl(seq, is.null))
}

# Remove rogues ----

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
  
  # Create OTU column for naming sequences as taxon-accession-gene
  metadata_with_seqs <- 
  metadata_with_seqs %>%
    # Make sure there are no spaces in taxon, accession, or gene
    verify(all(str_detect(taxon, " ", negate = TRUE))) %>%
    verify(all(str_detect(accession, " ", negate = TRUE))) %>%
    verify(all(str_detect(gene, " ", negate = TRUE))) %>%
    mutate(otu = glue("{taxon}-{accession}-{gene}"))
  
  # Extract sequences from metadata and rename
  pterido_seqs <- do.call(c, metadata_with_seqs$seq)
  names(pterido_seqs) <- metadata_with_seqs$otu
  
  # Make sure that went OK
  assertthat::assert_that(is.list(pterido_seqs))
  assertthat::assert_that(inherits(pterido_seqs, "DNAbin"))
  assertthat::assert_that(all(names(pterido_seqs) == metadata_with_seqs$otu))
  
  # Blast to exclude rogues
  
  # Make temporary working dir for BLAST functions.
  blast_dir <- fs::path(
      tempdir(), 
      # Use hash as unique name for folder
      digest::digest(pterido_seqs)
  )

  if(fs::dir_exists(blast_dir)) fs::dir_delete(blast_dir)

  fs::dir_create(blast_dir)
  
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
  if(fs::dir_exists(blast_dir)) fs::dir_delete(blast_dir)
  
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
  
  ### Check for monotypic families within each gene ###
  # These by default will match to a different (non-self)
  # family, so are not valid for checking rogues.
  exclude_from_rogues <-
    metadata_with_seqs %>% 
    dplyr::mutate(genus = stringr::str_split(taxon, "_") %>% purrr::map_chr(1)) %>%
    dplyr::left_join(
      select(ppgi, genus, family), 
      by = "genus") %>%
    assert(not_na, genus, family) %>%
    dplyr::add_count(family, gene) %>%
    dplyr::filter(n == 1) %>%
    dplyr::mutate(
      otu = glue("{taxon}-{accession}-{gene}") %>% stringr::str_replace_all(" ", "_")
    ) %>%
    select(otu, family, gene)
  
  # Group small families by order
  # to avoid false-positives
  ppgi <- mutate(ppgi, family = case_when(
    order == "Cyatheales" ~ "Cyatheales",
    suborder == "Saccolomatineae" ~ "Saccolomatineae",
    TRUE ~ family
  ))
  
  # Make list of rogue sequences (accessions) to exclude
  # Start with blast results
  blast_results %>%
    # Exclude monotypic families
    dplyr::anti_join(exclude_from_rogues, by = c(qseqid = "otu")) %>%
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
    # Add information for accession and gene (make sure OTU has expected number of dashes first)
    verify(all(str_count(qseqid, "-") == 2)) %>%
    mutate(
      accession = str_split(qseqid, "-") %>% map_chr(2),
      gene = str_split(qseqid, "-") %>% map_chr(3)) 
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

#' Select a set of GenBank genes by joining on common vouchers
#' 
#' This assumes genes comprise: "rbcL", "atpA", "atpB", "rps4". Filters
#' to only one set of sequences per taxon (species), prioritizing in order:
#' - those with rbcL + other genes
#' - those with rbcL
#' - those with the greatest combined length of other genes
#' 
#' @param genbank_seqs_tibble Tibble; metadata for GenBank sequences
#' including columns for `gene`, `species`, `specimen_voucher`, 
#' `accession` (GenBank accession), and `length` (gene length in bp)
#'
#' @return Tibble in wide format joining genes based on species + voucher
#' 
select_genbank_genes <- function (genbank_seqs_tibble) {
  
  ### Some pre-processing ### 
  # extract voucher specimen, count number of bases in sequence
  genbank_seqs_tibble_clean <- 
    genbank_seqs_tibble %>%
    # Filter to only records with specimen voucher
    select(subtype, subname) %>%
    unique() %>%
    filter(str_detect(subtype, "specimen_voucher")) %>%
    mutate(
      # Metadata are contained in two columns
      # `subtype` is column name of metadata separated by `|`
      # `subname` is content of metadata separated by `|`
      # `subtype` may not be unique for each record (could have multiple vouchers)
      # here, extract out `specimen_voucher`
      meta_names = str_split(subtype, "\\|"),
      meta_data = str_split(subname, "\\|"),
      meta_sel_index = map(meta_names, ~which(. == "specimen_voucher")),
    ) %>%
    unnest(meta_sel_index) %>%
    mutate(
      specimen_voucher = map2(meta_data, meta_sel_index, ~magrittr::extract(.x, .y))
    ) %>%
    unnest(specimen_voucher) %>%
    select(subtype, subname, specimen_voucher) %>%
    right_join(genbank_seqs_tibble, by = c("subtype", "subname")) %>%
    # Extract actual length of sequence (different from `slen` in gb metadata)
    mutate(length = map_dbl(seq, ~length(.[[1]]))) %>%
    select(gene, taxon, specimen_voucher, accession, length) %>%
    # In some cases, there are multiple vouchers including `s.n.` and 
    # a numbered voucher. Use only the numbered voucher.
    add_count(gene, taxon, accession) %>%
    filter(!(n > 1 & str_detect(specimen_voucher, "s\\.n\\."))) %>%
    # Now there should be a unique combination of `gene`, `taxon`, `accession` per row
    assert_rows(col_concat, is_uniq, gene, taxon, accession) %>%
    select(-n)
  
  ### Convert to wide format ###
  gb_data <-
    genbank_seqs_tibble_clean %>%
    group_by(gene, specimen_voucher, taxon) %>%
    arrange(desc(length)) %>%
    slice(1) %>%
    ungroup() %>%
    # Convert to wide format, joining on voucher
    # - first split into a list of dataframes by gene
    group_by(gene) %>%
    group_split() %>%
    # For each dataframe, convert to wide format and
    # rename "length" and "accession" columns by gene
    map(
      ~pivot_wider(.,
                   id_cols = c("taxon", "specimen_voucher"),
                   names_from = gene,
                   values_from = c("length", "accession")
      )
    ) %>%
    # Join the gene sequences by taxon + voucher
    reduce(full_join, by = c("taxon", "specimen_voucher"), na_matches = "never") %>%
    mutate_at(vars(contains("length")), ~replace_na(., 0)) %>%
    # Add total length of all genes (need to sum row-wise)
    # https://stackoverflow.com/questions/31193101/how-to-do-rowwise-summation-over-selected-columns-using-column-index-with-dplyr
    mutate(total_length = pmap_dbl(select(., contains("length")), sum))
  
  ### Select final sequences ###
  
  # Highest priority taxon: those with rbcL and at least one other gene.
  # Choose best specimen per taxon by total sequence length
  rbcl_and_at_least_one_other <-
    gb_data %>% 
    filter(!is.na(accession_rbcL)) %>%
    filter_at(
      vars(accession_atpA, accession_atpB, accession_rps4), 
      any_vars(!is.na(.))) %>%
    group_by(taxon) %>%
    arrange(desc(total_length)) %>%
    slice(1) %>%
    ungroup
  
  # Next priority: anything with rbcL only
  rbcl_only <-
    gb_data %>%
    filter(!taxon %in% rbcl_and_at_least_one_other$taxon) %>%
    filter(!is.na(accession_rbcL)) %>%
    group_by(taxon) %>%
    arrange(desc(total_length)) %>%
    slice(1) %>%
    ungroup
  
  # Next priority: any other taxon on basis of total seq. length
  other_genes <-
    gb_data %>%
    filter(!taxon %in% rbcl_and_at_least_one_other$taxon) %>%
    filter(!taxon %in% rbcl_only$taxon) %>%
    group_by(taxon) %>%
    arrange(desc(total_length)) %>%
    slice(1) %>%
    ungroup
  
  # Combine into final list
  bind_rows(
    rbcl_and_at_least_one_other,
    rbcl_only,
    other_genes
  ) %>%
    assert(is_uniq, taxon)
}

#' Filter out species in plastome data from Sanger data
#'
#' @param plastome_genes_unaligned List of unaligned plastome genes
#' @param plastome_metadata_renamed Plastome metadata after standardizing names
#' @param sanger_accessions_selection Selected GenBank (Sanger) sequences to use
#' @param filter Logical; should filtering be done or not?
#'
#' @return sanger_accessions_selection, with species in plastome data filtered out (if `filter` = TRUE) 
#' 
filter_out_plastome_species <- function (plastome_genes_unaligned, plastome_metadata_renamed, sanger_accessions_selection, filter) {
  
  # Don't do any filtering if `filter` is false
  if(filter == FALSE) return (sanger_accessions_selection)
  
  ### Make list of species in selected plastome genes ###
  # (not the same as plastome_metadata_renamed, since that was ALL plastomes, and
  # we need a list of species names of only the plastomes that will actually be used)
  plastid_genes_unaligned_species_list <-
    map_df(plastome_genes_unaligned, ~names(.) %>% tibble(accession = .), .id = "gene") %>%
    left_join(plastome_metadata_renamed, by = "accession") %>%
    select(gene, accession, species) %>%
    assert(not_na, species)
  
  # Filter out species from GenBank (Sanger) accessions that are already in plastome data
  anti_join(
    sanger_accessions_selection, 
    plastid_genes_unaligned_species_list,
    by = "species"
  )
  
}

# Download plastomes from GenBank ----

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
  
  # get the results
  extracted_genes <- map(genes, ~extract_sequence(gb_entry, .)) %>%
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
#' and genes absent from > 50% of sequences used in select_plastome_seqs()
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
#' and genes absent from > 50% of sequences used in select_plastome_seqs()
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
select_plastome_seqs <- function (plastid_seq_list, plastome_metadata, filter_by = c("species", "genus", "voucher")) {
  
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
#' @param plastome_selection Tibble of selected plastomes to use.
#'
#' @return List of unaligned genes.
#' 
extract_seqs_by_gene <- function (plastid_seq_list, plastome_selection) {
  plastid_seq_list %>%
    magrittr::extract(unique(plastome_selection$plastid_seq_name)) %>%
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

#' Write out a nexus block of gene start and end positions for IQTREE
#' 
#' Output nexus file example: http://www.iqtree.org/workshop/data/turtle.nex
#'
#' @param gene_list List of genes in the aligment
#' @param nexus_file Path to write out nexus file
#'
#' @return A vector of gene start and end positions
#' 
write_nexus_gene_block <- function (gene_list, nexus_file) {
  
  gene_locs <-
    tibble::tibble(
      gene = names(gene_list),
      length = purrr::map_dbl(gene_list, ncol)
    ) %>%
    dplyr::mutate(
      end = cumsum(length),
      start = end - length + 1) %>%
    dplyr::mutate(
      nexus_line = glue::glue("  charset {gene} = {start}-{end};")
    ) %>%
    dplyr::pull(nexus_line)
  
  cat(
    "#nexus", "begin sets;", gene_locs, "end;",
    file = nexus_file, sep = "\n"
  )
  
  gene_locs
  
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
      query == "Psilotum nudum" ~ "Psilotum nudum (L.) P. Beauv.",
      query == "Botrychium sp. ternatum/japonicum" ~ "Sceptridium ternatum (Thunb.) Lyon",
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
    select(query, resolved_name = species, scientificName)
  
  plastome_metadata %>%
    left_join(select(name_resolution_results_fixed, species = query, resolved_name, scientificName), by = "species") %>%
    mutate(species = resolved_name) %>%
    assert(not_na, species) %>%
    assert(not_na, scientificName) %>%
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
#' @param plastome_genes_unaligned Unaligned plastid genes from plastomes
#' @param plastome_metadata_renamed Plastome metadata with resolved names
#' @param pterido_rbcl_clean_seqs Unaligned rbcL sequences from genbank with
#' resolved names
#'
#' @return List; first item is the DNA sequences (list of class DNAbin); second
#' item is a dataframe with GenBank accession numbers for each species.
#' @export
#'
#' @examples
combine_rbcL <- function(plastome_genes_unaligned, plastome_metadata_renamed, pterido_rbcl_clean_seqs) {
  
  # Extract rbcL from the plastome genes list, and rename by species
  plastome_rbcL <- plastome_genes_unaligned[["rbcL"]]
  
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

#' Join parts of a DNA sequence together
#' using the range description in a GenBank flatfile
#' 
#' Helper function used by fetch_genes_from_plastome()
#'
#' @param seq Character vector of length 1; complete DNA sequence in format, e.g.,
#' "atgca"
#' @param ranges Character vector of length 1; ranges of the DNA sequence to join
#' together. Formatted as e.g., "1..3,8..10,20..40". No spaces allowed.
#'
#' @return Character vector of length 1
#' 
join_cds <- function(seq, ranges) {
  
  # Make sure only allowed characters are present
  assertthat::assert_that(
    str_detect(ranges, "[^0-9|..|,|<|>]") == FALSE,
    msg = "CDS range not properly specified"
  )
  
  # Convert range description string into list,
  # one item per range
  ranges_list <- 
    # Remove any characters other than '..', ',' or 0-9
    # (such as '<' or '>')
    str_remove_all(ranges, "[^0-9|..|,]") %>%
    str_split(",") %>% 
    magrittr::extract2(1) %>% 
    map(~str_split(., "\\.\\.") %>% magrittr::extract2(1))
  
  # Extract starts and ends of each range
  starts <- map_chr(ranges_list, 1)
  ends <-  map_chr(ranges_list, 2)
  
  # Convert sequence into character vector with one character each
  seq_vec <- str_split(seq, "")[[1]]
  
  # Extract the selected ranges, squash back together
  map2(starts, ends, ~magrittr::extract(seq_vec, .x:.y)) %>%
    unlist() %>%
    str_flatten()
}

#' Fetch a set of coding genes from (typically) a plastome
#'
#' @param accession GenBank accession number of the plastome
#' @param target_genes Character vector: list of target genes to extract
#' @param limit_missing Maximum number of missing amino acids to allow in the
#' coding sequence; those exceeding this number will be dropped
#'
#' @return List including
#' - dna: The DNA coding sequence of the gene, list of class "DNAbin"
#' - aa: The amino acid sequence of the gene, list of class "AAbin"
#' - duplicates: Dataframe of genes that appear more than once in the
#' plastome and are excluded from 'dna' and 'aa'.
#'
fetch_genes_from_plastome <- function (accession, target_genes, limit_missing = 10) {
  
  # Download and read in GenBank flatfile ---
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
  
  # Parse GenBank record
  # suppress "Sample 1 in 1 done" message
  read_gb <- quietly(read.gb::read.gb)
  gb_parsed <- read_gb(temp_file) %>% 
    magrittr::extract2("result") %>%
    flatten()
  
  fs::file_delete(temp_file)
  
  # Extract full DNA sequence ---
  dna <- gb_parsed[["ORIGIN"]]
  
  # Make sure it only includes IUPAC bases
  unique_bases <- dna %>% str_split("") %>% magrittr::extract2(1) %>% unique()
  
  iupac <- c("A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N") %>%
    c(., str_to_lower(.))
  
  assertthat::assert_that(
    all(unique_bases %in% iupac),
    msg = "DNA sequence includes non-IUPAC bases")
  
  # Extract CDS with DNA sequence and translations ---
  cds_with_dups <- 
    # Extract CDS from parsed flat-file (includes translation and DNA range)
    gb_parsed[["FEATURES"]][names(gb_parsed[["FEATURES"]]) == "CDS"] %>%
    map_df(~pivot_wider(., values_from = Qualifier, names_from = Location)) %>%
    # Filter to only target genes
    filter(gene %in% target_genes) %>% 
    select(gene, CDS, translation) %>%
    group_by(gene) %>%
    # Count number of entries per gene for each CDS
    add_count() %>%
    ungroup # Be sure to ungroup before using assertr!
  
  duplicates <- filter(cds_with_dups, n > 1)
  
  cds_with_all_aa <-
    # Filter out any genes with >1 CDS
    cds_with_dups %>%
    filter(n == 1) %>%
    select(-n) %>%
    # Verify that parentheses are paired, and there are no more than two nested parentheses
    mutate(
      n_open_paren = str_count(CDS, "\\("),
      n_close_paren = str_count(CDS, "\\(")) %>%
    verify(n_open_paren == n_close_paren) %>%
    verify(n_open_paren <= 2) %>%
    select(-n_open_paren, -n_close_paren) %>%
    # Get rid of "join" and closing parentheses
    mutate(CDS = str_remove_all(CDS, "join\\(|\\)")) %>%
    # Parse CDNA range
    separate(CDS, into = c("complement", "range"), sep = "\\(", fill = "left") %>%
    mutate(range = str_remove_all(range, ",$")) %>%
    # Get rid of extra spaces and commas in translation
    mutate(translation = str_remove_all(translation, " |,")) %>%
    # Filter remaining by number of missing residues
    mutate(
      aa_len = nchar(translation),
      # Set sequence name as accession-gene (format required by hybpiper)
      seq_name = paste(all_of(accession), gene, sep = "-"),
      n_aa_missing = str_count(translation, "X")
    ) %>%
    filter(n_aa_missing < limit_missing) %>%
    # Add DNA sequence
    mutate(dna_seq = map_chr(range, ~join_cds(all_of(dna), .x))) %>%
    # Check that DNA sequence length matches AA length
    # (should be within 1 AA of 1 codon)
    mutate(
      dna_len = nchar(dna_seq),
      expect_aa_low = floor((dna_len - 3) / 3),
      check_aa_low = aa_len >= expect_aa_low,
      expect_aa_high = ceiling((dna_len + 3) / 3),
      check_aa_high = aa_len <= expect_aa_high
    )
  
  cds <- cds_with_all_aa %>%
    filter(check_aa_low, check_aa_high)
  
  aa_length_errors <- anti_join(cds_with_all_aa, cds, by = "gene")
  
  # Convert AA and DNA to APE format
  aa <- ape::as.AAbin(as.list(cds$translation)) %>%
    set_names(cds$seq_name)
  
  dna <-
    cds$dna_seq %>%
    map(~str_split(., "")) %>%
    map(ape::as.DNAbin) %>%
    set_names(cds$seq_name)
  
  # Reverse-complement those DNA sequences that need it
  dna[which(cds$complement == "complement")] <- dna[which(cds$complement == "complement")] %>% map(ape::complement)
  
  list(
    dna = dna,
    aa = aa,
    duplicates = duplicates,
    length_errors = aa_length_errors
  )
  
}


#' Combine sequences for target genes from GenBank with sequences
#' extracted from plastomes
#'
#' @param raw_fasta_all_genes Tibble of DNA sequence data downloaded from GenBank
#' @param sanger_accessions_selection Final selection of accessions to use after removing rogues
#' @param plastome_genes_unaligned List of unaligned genes extracted from plastomes
#'
#' @return List of class DNAbin
#' 
combine_sanger_with_plastome <- function (
  raw_fasta_all_genes,
  sanger_accessions_selection,
  plastome_genes_unaligned
) {
  
  # Convert final selected GenBank (Sanger) accessions to long format
  final_gb_accessions <-
    sanger_accessions_selection %>%
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
  common_gene_names <- intersect(names(genbank_genes_unaligned), names(plastome_genes_unaligned)) 
  
  # - Make a list of accessions for genes in common between Sanger and plastome sequences
  common_genes <-
    common_gene_names %>%
    map(~c(genbank_genes_unaligned[[.]], plastome_genes_unaligned[[.]])) %>%
    set_names(common_gene_names)
  
  # - Make a list of accessions in Sanger genes only (likely zero, but for completeness' sake)
  genbank_genes_only <- genbank_genes_unaligned %>%
    magrittr::extract(setdiff(names(genbank_genes_unaligned), names(plastome_genes_unaligned)))
  
  # - Make a list of accessions in plastome sequences only
  plastome_genes_only <- plastome_genes_unaligned %>%
    magrittr::extract(setdiff(names(plastome_genes_unaligned), names(genbank_genes_unaligned)))
  
  # Combine the lists
  c(common_genes, genbank_genes_only, plastome_genes_only)
  
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

# fastp ----

#' Parse a JSON file output by fastp into a dataframe
#'
#' @param file Path to the JSON file
#' @param sample_name Sample name (will be automatically
#' detected from  files named like 's_1_1_sequence_trim_report.json'
#' if not provided)
#'
#' @return Dataframe (tibble)
#' 
parse_fastp_sum <- function (file, sample_name = NULL) {
  
  # Parse JSON with jsonlite package
  data <- jsonlite::fromJSON(file)
  
  if(is.null(sample_name)) sample_name <- fs::path_file(file) %>% str_remove_all("_sequence_trim_report.json")
  
  # Extract data of interest as tibbles
  before_filtering <- as_tibble(data[["summary"]][["before_filtering"]]) %>%
    rename_with(~paste0(.x, "_before"))
  
  after_filtering <- as_tibble(data[["summary"]][["after_filtering"]]) %>%
    rename_with(~paste0(.x, "_after"))
  
  filtering_result <- as_tibble(data[["filtering_result"]])
  
  duplication_rate <- tibble(duplication_rate = data[["duplication"]][["rate"]])
  
  # Combine the tibbles
  bind_cols(
    sample = sample_name,
    before_filtering,
    after_filtering,
    filtering_result,
    duplication_rate
  )
  
}

#' Run fastp
#' 
#' Trims adapters and low-quality bases on fastp default settings, assuming
#' paired-end input.
#' 
#' For more about fastp, see:
#' https://github.com/OpenGene/fastp
#' 
#' Trimmed reads along with an html report will be be written to `out_dir`.
#'
#' @param sample String indicating unique sample name (other than suffix specifying
#' forward and reverse reads)
#' @param data_dir Location of samples (actual samples may be nested within this)
#' @param out_dir Location to write out filtered reads
#' @param f_suffix Suffix specifying forward reads
#' @param r_suffix Suffix specifying reverse reads
#'
#' @return Tibble; stats from trimming reads
#' @export
#'
#' @examples
fastp <- function(sample, data_dir = "data_raw", out_dir = "intermediates/fastp/", f_suffix = "R1_001", r_suffix = "R2_001") {
  
  # Find paths to raw F and R reads
  raw_forward <- list.files(data_dir, full.names = TRUE, pattern = glue(".*{sample}.*{f_suffix}"), recursive = TRUE)
  raw_reverse <- list.files(data_dir, full.names = TRUE, pattern = glue(".*{sample}.*{r_suffix}"), recursive = TRUE)
  
  # Verify that input paths are OK
  assertthat::assert_that(length(raw_forward) == 1, msg = "Sample does not match exactly one forward raw read")
  assertthat::assert_that(length(raw_reverse) == 1, msg = "Sample does not match exactly one reverse raw read")
  assertthat::assert_that(assertthat::is.readable(raw_forward))
  assertthat::assert_that(assertthat::is.readable(raw_reverse))
  
  # Specify temporary file for writing trimming report in json format
  temp_json <- fs::path(tempdir(), glue("{sample}.json"))
  
  # Set up fastp arguments
  args <- c(
    "-i", raw_forward,
    "-I", raw_reverse,
    "-o", glue::glue("{out_dir}/{sample}_R1.fastq"),
    "-O", glue::glue("{out_dir}/{sample}_R2.fastq"),
    "-j", temp_json,
    "-h", glue::glue("{out_dir}/{sample}.html")
  )
  
  # Run fastp
  processx::run("fastp", args)
  
  # Parse json report into tibble
  res <- parse_fastp_sum(temp_json, sample)
  
  # Cleanup
  fs::file_delete(temp_json)
  
  # Return tibble
  res
  
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
  Sys.sleep(1)
  
  assertthat::assert_that(assertthat::is.string(accession))
  assertthat::assert_that(assertthat::is.string(rank))
  
  # Check for ENTREZ KEY
  assertthat::validate_that(
    nchar(Sys.getenv("ENTREZ_KEY")) > 0,
    msg = "Warning: ENTREZ_KEY not detected. See taxize::key_helpers()."
  )
  
  taxonomy <-
    taxize::genbank2uid(id = accession) %>%
    taxize::classification(db = "ncbi", max_tries = 10) %>%
    purrr::flatten()
  
  if(rank %in% taxonomy$rank) return(taxonomy$name[taxonomy$rank == rank])
  
  NA_character_
  
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
    assert(not_na, species) %>%
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
  
  # Specify paths in new environment
  new_env <- Sys.getenv()
  new_env["PATH"] <- paste(Sys.getenv("PATH"), "/apps/SPAdes/3.13.0/bin/:/apps/HybPiper/:/apps/partitionfinder/2.1.1/", sep = ":")
  
  # Run command
  processx::run(
    "reads_first.py", 
    hybpiper_arguments, 
    wd = wd, echo = echo, env = new_env)
}

# Retrieve read fragments from hybpiper ----

#' Extract the consensus sequence of short reads sorted by HybPiper
#' 
#' Requires bbmap and kindel to be installed and on the PATH.
#'
#' @param sample Name of sample. Should be the same as the name of the sample
#' used as input to Hypiper. Assumes hybiper has been run, with output in the
#' "intermediates/hybpiper" folder.
#' @param plastid_targets List of lists, one per target gene. Each list includes
#' - dna: The DNA coding sequence of the gene, list of class "DNAbin"
#' - aa: The amino acid sequence of the gene, list of class "AAbin"
#' - duplicates: Dataframe of genes that appear more than once in the
#' plastome and are excluded from 'dna' and 'aa'.
#' Output of map(accessions, ~fetch_genes_from_plastome(., wei_genes))
#' @param ... Other arguments not used by this function but meant for tracking
#' with `drake`.
#'
#' @examples
#' source("_skimming_drake.R")
#' # Load plastid_targets
#' loadd(plastid_targets, cache = skimming_cache)
#' # sample with fewest reads
#' sample <- "UFG_393201_P03_WC08"
#' get_hybpip_consensus("UFG_393201_P03_WC08", plastid_targets)
get_hybpip_consensus <- function (sample, plastid_targets, ...) {
  
  # Make temp working folder
  temp_dir <- fs::dir_create(fs::path(tempdir(), digest::digest(sample)))
  
  # Load list of all DNA targets, format
  # as tibble named by accession-gene
  target_dna <-
    tibble(dna = transpose(plastid_targets)[["dna"]] %>% flatten) %>%
    mutate(
      name = names(dna),
      gene = str_split(name, "-") %>% map_chr(2),
      dna = map2(dna, name, ~set_names(.x, .y))
    )
  
  # Fetch the names of target sequences matched to the sample by hybpiper
  ref_names <- list.files(
    paste0("intermediates/hybpiper/", sample), 
    pattern = "baits.fasta", 
    full.names = TRUE, recursive = TRUE) %>%
    map(ape::read.FASTA) %>%
    map_chr(names)
  
  # Construct pseudo reference genome based on the matched targets
  pseudo_ref_genome <- filter(target_dna, name %in% ref_names) %>%
    pull(dna) %>%
    set_names(nm = NULL) %>%
    jntools::flatten_DNA_list()
  
  # Write out pseudo reference genome
  ape::write.FASTA(pseudo_ref_genome, fs::path(temp_dir, "ref.fasta"))
  
  # Load all the sorted reads from HybPiper for this sample, combine
  short_reads <- list.files(
    paste0("intermediates/hybpiper/", sample), 
    pattern = "interleaved.fasta", 
    full.names = TRUE, recursive = TRUE) %>%
    map(ape::read.FASTA) %>%
    jntools::flatten_DNA_list()
  
  # Write out sorted reads for mapping
  ape::write.FASTA(short_reads, fs::path(temp_dir, "reads.fasta"))
  
  # Run bbmap: maps reads to reference
  args <- c(
    # input short reads
    glue('in={fs::path(temp_dir, "reads.fasta")}'),
    # reference genome to map reads
    glue('ref={fs::path(temp_dir, "ref.fasta")}'),
    # output path to write aligned reads (SAM)
    glue('out={fs::path(temp_dir, "align.sam")}'),
    # use one thread
    "t=1",
    # reserve 1gb memory
    "-Xmx1g",
    # don't write out the index file to disk
    "-nodisk"
  )
  
  bbmap_res <- processx::run("bbmap", args)
  
  # Process bbmap stderr to get table of reads mapped, 
  # number of input reads, number of input bases
  bbmap_stderr_raw <- write_lines(bbmap_res$stderr, fs::path(temp_dir, "bbmap.stderr"))
  bbmap_stderr_lines <- read_lines(bbmap_stderr_raw)
  skip_val <- str_detect(bbmap_stderr_lines, "Read 1 data") %>% which %>% magrittr::subtract(1)
  max_val <- str_detect(bbmap_stderr_lines, "Total time") %>% which %>% magrittr::subtract(skip_val + 6)
  
  map_stats <- read_tsv(bbmap_stderr_raw, skip = skip_val, n_max = max_val) %>%
    janitor::clean_names() %>%
    rename(category = read_1_data) %>%
    mutate(category = str_remove_all(category, ":"))
  
  n_input_reads <- str_match_all(bbmap_res$stderr, "Reads Used: +\t([0-9]+)") %>% purrr::pluck(1,2)
  n_input_bp <- str_match_all(bbmap_res$stderr, "Reads Used:[^\\(]+\\(([0-9]+) ") %>% purrr::pluck(1,2)
  
  # Run kindel: extracts consensus from alignment
  args <- c(
    # kindel subcommand: produce consensus
    "consensus",
    # input alignment file
    fs::path(temp_dir, "align.sam"),
    # set required read depth to 1 (the minimum)
    "--min-depth", "1",
    # trim ambiguous bases from ends
    "-t"
  )
  
  kindel_res <- processx::run("kindel", args)
  
  # Write out standard output so it can be read in with ape::read.FASTA
  con_raw <- kindel_res$stdout %>%
    # remove some text from the sequence names added by kindel
    str_remove_all("_cns <unknown description>")
  
  write_lines(con_raw, fs::path(temp_dir, "con.fasta"))
  
  # Read back in with ape
  consensus <- ape::read.FASTA(fs::path(temp_dir, "con.fasta"))
  
  # Reformat sequence names to "sample-gene"
  new_names <- paste(
    sample,
    names(consensus) %>% str_match("-(.+)$") %>% magrittr::extract(,2),
    sep = "-"
  )
  
  names(consensus) <- new_names
  
  # Combine results into tibble
  results <- tibble(
    consensus = list(consensus),
    map_stats = list(map_stats),
    n_input_reads = n_input_reads,
    n_input_bp = n_input_bp,
    sample = sample
  )
  
  # Delete temporary files
  fs::dir_delete(temp_dir)
  
  return(results)
  
}

# Reporting ----

#' Make a table mapping accession numbers to species and genes
#'
#' @param plastid_genes_aligned_trimmed List of plastid sequences used for phylogenetic
#' analysis. Aligned and trimmed, named by GenBank accession.
#' @param sanger_seqs_names_resolved Metadata with resolved species names and accession
#' numbers for Sanger sequences, including some not used in phylo. analysis.
#' @param plastome_metadata_renamed Metadata with resolved species names and accession
#' numbers for plastome sequences, including some not used in phylo. analysis.
#' @param target_genes Character vector of Sanger genes included in phylo. analysis
#'
#' @return Tibble
#' 
make_acc_ref_table <- function(
  plastid_genes_aligned_trimmed,
  sanger_seqs_names_resolved,
  plastome_metadata_renamed,
  target_genes
) {
  
  # Get a vector of all accession numbers in the final plastid alignment
  plastid_accs <- map(plastid_genes_aligned_trimmed, rownames) %>%
    set_names(NULL) %>%
    unlist()
  
  # combine sanger and plastome metata data
  sanger_plastome_dat <-
    bind_rows(
      select(sanger_seqs_names_resolved, accession, species, sci_name = scientificName, gene),
      transmute(plastome_metadata_renamed, accession, species, sci_name = scientificName, gene = "full_plastome")
    )
  
  # Build table of accession numbers, species, gene, and "title" for full plastome sequences
  # - start with accession numbers used for phy. analysis.
  tibble(accession = plastid_accs) %>%
    mutate(accession = str_remove_all(accession, "_R_")) %>%
    unique() %>%
    left_join(sanger_plastome_dat, by = "accession") %>%
    # Make sure we aren't missing anything
    assert(not_na, species, sci_name, gene) %>%
    # Make sure the combination of gene/accession is unique
    assert_rows(col_concat, is_uniq, accession, gene) %>%
    # Convert to wide form
    pivot_wider(names_from = "gene", values_from = "accession", species:sci_name) %>%
    arrange(species) %>%
    # Rearrange columns
    select(species, sci_name, all_of(target_genes), full_plastome)
}

#' Make a table of gene partitions
#'
#' @param gene_list List of gene alignments that are
#' used for concatenation
#'
#' @return Tibble of start and end positions of each gene
#' in the concatenated data
#' 
make_gene_part_table <- function(gene_list) {
  tibble(
    gene = names(gene_list),
    length = map_dbl(gene_list, ncol)) %>%
    mutate(
      end = cumsum(length),
      start = end - (length - 1)
    ) %>%
    select(gene, start, end) %>%
    assert(is_uniq, gene)
}

# Managing data ----

#' Make a zipped archive of raw data
#'
#' @param version Version number to assign to the archive
#' @param metadata Tibble with the following metadata:
#' - target: name of target for raw data file in drake plan (only if that file
#' is read in individually as a target)
#' - file: file name
#' - path: path to file (relative to this project)
#' - copy: Boolean; should the file by copied into the archive?
#' - hash: MD5 hash of the file
#' @param out_path Directory to save the zip archive
#'
#' @return Status of the `zip` program after running; externally, the data
#' will be saved to `<version>.zip` in `out_path`.
#'
archive_raw_data <- function (version, metadata, out_path) {
  
  # Specify location of zip archive
  archive <- paste0(fs::path(out_path, version), ".zip")
  
  # Make sure the zip archize doesn't already exist
  assertthat::assert_that(
    !fs::file_exists(archive), 
    msg = glue::glue("{archive} already exists")
  )
  
  # Create a temp dir for writing out csv file with metadata and README
  temp_dir <- tempdir()
  
  # Write out csv of metadata (MD5 checksums)
  readr::write_csv(metadata, fs::path(temp_dir, "md5_checksums.csv"))
  
  # Make vector of raw data files to copy into zip folder
  raw_data_to_copy <- metadata %>%
    filter(copy == TRUE) %>%
    # When the files get placed into the zip folder, they won't
    # be in folders, so make sure each has a unique name.
    assertr::assert(assertr::is_uniq, file) %>%
    pull(path) %>%
    c(fs::path(temp_dir, "md5_checksums.csv"))
  
  # zip the files
  zip_results <- utils::zip(
    zipfile = archive, 
    files = raw_data_to_copy,
    flags = "-r9Xj")
  
  # Cleanup
  fs::file_delete(fs::path(temp_dir, "md5_checksums.csv"))
  
  # Returnc
  zip_results
  
}

# Taxonomic name resolution ----

#' Make a dataframe with taxonomic names
#'
#' @param taxa Character vector; taxon names to be parsed by taxon-tools `parsenames`.
#' Missing values not allowed. Must all be unique.
#'
#' @return Dataframe with two columns: `id` and `name`
#' 
#' @examples
#' tt_make_name_df("Foogenus x barspecies var. foosubsp (L.) F. Bar")
tt_make_name_df <- function(taxa) {
  
  assertthat::assert_that(all(assertr::is_uniq(taxa)))
  assertthat::assert_that(assertthat::noNA(taxa))
  
  # Format input names as data frame with unique ID
  # ID is combination of first 4 chars of hash of the input (taxa), followed by "-" and integer
  taxa_tbl <- data.frame(name = taxa)
  taxa_tbl$id <- 1:nrow(taxa_tbl)
  taxa_tbl$id <- paste(substr(digest::digest(taxa), 1, 4), taxa_tbl$id, sep = "-")
  
  taxa_tbl[, c("id", "name")]
}

#' Parse taxonomic names with (and for) taxon-tools
#' 
#' Requires [taxon-tools](https://github.com/camwebb/taxon-tools) to be installed.
#' 
#' Parses scientific names into their component parts (genus, species, variety, author, etc).
#'
#' @param taxa Character vector; taxon names to be parsed by taxon-tools `parsenames`.
#' Missing values not allowed. Must all be unique.
#' @param file Character vector of length 1; path to write parsed names (optional).
#'
#' @return A dataframe including the following columns.
#' - id: A unique ID number assigned to the input name
#' - name: The input name
#' - genus_hybrid_sign: Hybrid sign for genus
#' - genus_name: Genus name
#' - species_hybrid_sign: Hybrid sign for species
#' - specific_epithet: Specific epithet (name)
#' - infraspecific_rank: Infraspecific rank
#' - infraspecific_epithet: Infraspecific epithet (name)
#' - author: Name of taxon
#' 
#' If `file` is specified, the parsed names will be written to that path in a format
#' for taxon-tools `matchnames`.
#' 
#' @examples 
#' tt_parse_names(c("Foogenus x barspecies var. foosubsp (L.) F. Bar", " "))
#'
tt_parse_names <- function(taxa) {
  
  # Check input: must be character vector, no NA values, all unique
  assertthat::assert_that(is.character(taxa))
  
  # Write out names formatted for parsing with taxon-tools to temp file
  # format: 
  # `id_num|taxon_name`
  # for example,
  # `x-234|Foogenus x barspecies var. foosubsp (L.) F. Bar`
  taxa_tbl <- tt_make_name_df(taxa)
  taxa_tbl$record <- paste(taxa_tbl$id, taxa_tbl$name, sep = "|")
  ref_taxa_txt_file <- tempfile(
    pattern = digest::digest(taxa),
    fileext = ".txt"
  )
  if(fs::file_exists(ref_taxa_txt_file)) fs::file_delete(ref_taxa_txt_file)
  writeLines(taxa_tbl$record, ref_taxa_txt_file)
  
  # Parse reference names with taxon tools
  ref_parsed <- processx::run("parsenames", ref_taxa_txt_file)
  if(fs::file_exists(ref_taxa_txt_file)) fs::file_delete(ref_taxa_txt_file)
  
  # Read in results of parsing, format as dataframe
  
  # The output is originally one record per line, with fields separated by '|' (pipe symbol)
  parsed_res <- data.frame(record = strsplit(ref_parsed[["stdout"]], "\n")[[1]])
  
  # Split these into separate columns
  name_parts <- c(
    "genus_hybrid_sign",
    "genus_name",
    "species_hybrid_sign",
    "specific_epithet",
    "infraspecific_rank",
    "infraspecific_epithet",
    "author"
  )
  
  parsed_res <- tidyr::separate(
    data = parsed_res, 
    col = record, 
    into = c("id", name_parts), 
    sep = "\\|", 
    fill = "right",
    remove = FALSE)
  
  # Fill in NA if that name part is missing
  parsed_res[parsed_res == ""] <- NA
  
  # Add "fail" column if all name parts are missing (couldn't be parsed properly)
  parsed_res$fail <- sapply(1:nrow(parsed_res), function(x) all(is.na(parsed_res[x, name_parts])))
  
  # Early exit if everything failed
  assertthat::assert_that(
    !all(parsed_res$fail == TRUE),
    msg = "No names could be successfully parsed")
  
  # Emit warning for failures
  if(sum(parsed_res$fail) > 0) {
    failed_ids <- parsed_res$id[parsed_res$fail == TRUE]
    failed_names <- paste(taxa_tbl$name[taxa_tbl$id %in% failed_ids], collapse = ", ") 
    warning(glue::glue("The following names could not be parsed and are excluded from results: {failed_names}"))
  }
  
  # Add back in original name
  parsed_res <- dplyr::left_join(
    parsed_res,
    dplyr::select(taxa_tbl, id, name), by = "id")
  
  # Remove failures, drop "fail" column
  parsed_res <- parsed_res[parsed_res$fail == FALSE, ]
  parsed_res$fail <- NULL
  
  # Return parsed names as dataframe
  parsed_res[, c("name", "id", name_parts)]
}

#' Load parsed taxon names
#'
#' @param path Path to taxon names that were parsed with taxon-tools
#' https://github.com/camwebb/taxon-tools
#'
#' @return Dataframe (tibble)
#' 
tt_load_names <- function(path) {
  suppressMessages(
    # Read in the parsed names as a dataframe
    readr::read_delim(
      path, delim = "|", 
      col_types = cols(.default = readr::col_character()),
      col_names = c(
        "id",
        "genus_hybrid_sign",
        "genus_name",
        "species_hybrid_sign",
        "specific_epithet",
        "infraspecific_rank",
        "infraspecific_epithet",
        "author"
      )
    )
  )
}

#' Write out parsed names for taxon-tools from a dataframe to a text file
#'
#' @param df Dataframe with parsed names
#' @param path Path to write dataframe
#'
#' @return Path to parsed names
tt_write_names <- function(df, path) {
  
  # Make vector of standard taxon-tools columns
  tt_col_names = c(
    "id",
    "genus_hybrid_sign",
    "genus_name",
    "species_hybrid_sign",
    "specific_epithet",
    "infraspecific_rank",
    "infraspecific_epithet",
    "author"
  )
  
  # Replace NA values with ""
  df <- dplyr::mutate(df, dplyr::across(dplyr::everything(), ~tidyr::replace_na(., "")))
  
  # Subset to only taxon-tools columns, in order
  df <- df[, tt_col_names]
  
  # taxon-tools uses pipe as separator
  df <- tidyr::unite(df, col = "text", dplyr::all_of(tt_col_names), sep = "|")
  
  # write out text
  writeLines(df$text, path)
  
  path
  
}

#' Match names using taxon-tools
#' 
#' Requires [taxon-tools](https://github.com/camwebb/taxon-tools) to be
#' installed.
#' 
#' `taxon-tools` matches names in two steps: first scientific names are
#' parsed into their component parts (genus, species, variety, author, etc). 
#' Next, names are fuzzily matched following rules using the component parts.
#' 
#' Parsing is fairly fast (much faster than matching) but can take some time if the
#' number of names is very large. If multiple queries will be made, it is recommended
#' to first parse the names to a file using `tt_parse_names()`, and use the file(s) as
#' input to `query` and/or `reference`.
#' 
#' One of either `reference` or `cache` must be provided.
#'
#' @param query Character vector; taxon names to be queried, or the path to a text
#' file with taxon names that have been parsed with `tt_parse_names()`.
#' Missing values not allowed. Must all be unique.
#' @param reference  Character vector; taxon names to use as reference, or the path to a text
#' file with taxon names that have been parsed with `tt_parse_names()`.
#' Missing values not allowed. Must all be unique.
#' @param max_dist Max Levenshtein distance to allow during fuzzy matching. 
#' (total insertions, deletions and substitutions)
#' @param match_no_auth Logical; If no author is given in the query and the name (without author) 
#' occurs only once in the reference, accept the name in the reference as a match. 
#' Default: to not allow such a match.
#' @param match_canon Logical; Allow a "canonical name" match if only the genus, species epithet, 
#' and infraspecific epithet (if present) match exactly. Default: to not allow such a match.
#' @param simple Logical; return the output in a simplified format with only the query
#' name, matched reference name, and match type.
#'
#' @return Dataframe
#' 
tt_match_names <- function(
  query, reference, 
  max_dist = 10, match_no_auth = FALSE, match_canon = FALSE, simple = FALSE
) {
  
  # Check input
  assertthat::assert_that(is.character(query) | inherits(query, "data.frame"))
  # assertthat::assert_that(assertthat::noNA(query))
  # assertthat::assert_that(all(assertr::is_uniq(query)))
  assertthat::assert_that(is.character(reference) | inherits(reference, "data.frame"))
  # assertthat::assert_that(assertthat::noNA(reference))
  # assertthat::assert_that(all(assertr::is_uniq(reference)))
  assertthat::assert_that(assertthat::is.number(max_dist))
  assertthat::assert_that(is.logical(match_no_auth))
  assertthat::assert_that(is.logical(match_canon))
  assertthat::assert_that(is.logical(simple))
  
  # Parse or load query names
  if(is.character(query)) {
    # Parse the names (adds 'name' column)
    query_parsed_df <- tt_parse_names(query)
  } else {
    # Or, names are already parsed
    query_parsed_df <- query
  }
  
  # Write out parsed names to temporary file
  query_parsed_txt <- tempfile(pattern = digest::digest(query), fileext = ".txt")
  if(fs::file_exists(query_parsed_txt)) fs::file_delete(query_parsed_txt)
  tt_write_names(query_parsed_df, query_parsed_txt)
  
  # Parse or load reference names
  if(is.character(reference)) {
    # Parse the names (adds 'name' column)
    ref_parsed_df <- tt_parse_names(reference)
  } else {
    # Or, names are already parsed
    ref_parsed_df <- reference
  }
  
  # Write out parsed names to temporary file
  ref_parsed_txt <- tempfile(pattern = digest::digest(reference), fileext = ".txt")
  if(fs::file_exists(ref_parsed_txt)) fs::file_delete(ref_parsed_txt)
  tt_write_names(ref_parsed_df, ref_parsed_txt)
  
  # Format argument flags
  if(match_no_auth) match_no_auth <- "-1" else match_no_auth <- NULL
  if(match_canon) match_canon <- "-c" else match_canon <- NULL
  
  # Specify temporary output file
  match_results_txt <- tempfile(pattern = digest::digest(c(query, reference)), fileext = ".txt")
  if(fs::file_exists(match_results_txt)) fs::file_delete(match_results_txt)
  
  # Run taxon-tools matchnames
  match_results <- processx::run(
    command = "matchnames", 
    args = c(
      "-a", query_parsed_txt,
      "-b", ref_parsed_txt,
      "-o", match_results_txt,
      "-e", max_dist,
      "-F", # no manual matching
      match_no_auth,
      match_canon
    )
  )
  
  # Read in results
  # Each line represents a single name from the query list (list A). 
  # Seventeen pipe-delimited (“|”) fields per row: 
  #  1. User ID code in list A, 
  #  2. Code in list B (if matched), 
  #  3. Match type (see codes below), 
  #  4-10. Parsed elements of name in list A. 
  #  11-17 (in same format as name input), Parsed elements of name in list B.
  matchnames_cols <- c(
    "id_query",
    "id_ref",
    "match_type",
    "genus_hybrid_sign_query",
    "genus_name_query",
    "species_hybrid_sign_query",
    "specific_epithet_query",
    "infraspecific_rank_query",
    "infraspecific_epithet_query",
    "author_query",
    "genus_hybrid_sign_ref",
    "genus_name_ref",
    "species_hybrid_sign_ref",
    "specific_epithet_ref",
    "infraspecific_rank_ref",
    "infraspecific_epithet_ref",
    "author_ref"
  )
  
  results <- data.frame(record = readLines(match_results_txt))
  
  results <- tidyr::separate(
    data = results, 
    col = record, 
    into = matchnames_cols, 
    sep = "\\|", 
    fill = "right",
    remove = TRUE)
  
  # Add back in the original search terms (query and reference)
  results <- dplyr::left_join(
    results,
    dplyr::select(query_parsed_df, id_query = id, query = name), 
    by = "id_query") 
  
  results <- dplyr::left_join(
    results,
    dplyr::select(ref_parsed_df, id_ref = id, reference = name), 
    by = "id_ref")
  
  if(simple == TRUE) results <- dplyr::select(results, query, reference, match_type)
  
  results
  
}

#' Fetch taxonomic data from the NCBI taxonomy database
#' 
#' This should only be run on a maximum of 200 IDs at a time,
#' or NCBI won't accept the query
#'
#' @param tax_ids Vector of taxonomy database IDs
#'
#' @return List (result of parsing XML output to a list)
#'
#' @examples
#' fetch_taxonomy_chunk("2726184")
fetch_taxonomy_chunk <- function(tax_ids) {
  rentrez::entrez_fetch(
    db = "taxonomy", 
    id = paste(tax_ids, collapse = "|"), 
    rettype="xml", 
    retmode = "text", 
    retmax = length(tax_ids)) %>%
    # Original output of rentrez is plain XML.
    # Parse this into a list, with one item for each of `tax_ids`
    stringr::str_split("\\n") %>%
    magrittr::extract2(1) %>%
    XML::xmlToList()
}

#' Parse a taxonomic data in a list downloaded from the NCBI taxonomy database
#'
#' @param record Single set of taxonomic data (one record) downloaded from
#' the NCBI taxonomy database with fetch_taxonomy_chunk()
#'
#' @return Dataframe with columns "taxid", "species", "scientific_name", and "accepted"
#' 
#' @examples
#' tax_data <- fetch_taxonomy_chunk("2726184")
#' parse_tax_list(tax_data[[1]])
parse_tax_list <- function(record) {
  
  # Each record has only 1 taxon id
  taxid <- record$TaxId
  # Each record has only 1 accepted species name 
  accepted_species <- record$ScientificName
  # "OtherNames" includes synonyms and scientific names in a nested list
  # Need to select by name directly because multiple items in the list may have the
  # same name (e.g., multiple "Synonym", etc)
  other_names <- record[names(record) == "OtherNames"] %>% purrr::flatten()
  # Confusingly, "synonyms" may sometimes contain the accepted name :/
  synonyms <- other_names[names(other_names) == "Synonym"] %>% unlist()
  sci_names <- other_names[names(other_names) == "Name"] %>% purrr::flatten()
  sci_names <- sci_names[names(sci_names) == "DispName"] %>% unlist()
  
  # The results should always have at least taxon ID and species
  species_dat <- tibble(taxid = taxid, species = accepted_species, accepted = TRUE)
  
  # Create empty tibbles to hold other name data
  acc_sci_names_dat <- tibble()
  syn_sci_names_dat <- tibble()
  
  # If other sci names are given, one of them should be the species name
  if(!is.null(sci_names)) 
    acc_sci_names_dat <- tibble(scientific_name = sci_names) %>%
    mutate(species = str_extract(scientific_name, accepted_species)) %>%
    filter(!is.na(species)) %>%
    mutate(accepted = TRUE)
  
  # If "synonym" and other sci names are given, one (or more) of them are the synonym
  if(!is.null(synonyms) & !is.null(sci_names)) 
    syn_sci_names_dat <- tibble(scientific_name = sci_names) %>%
    mutate(species = str_extract(scientific_name, paste(synonyms, collapse = "|"))) %>%
    filter(!is.na(species)) %>%
    mutate(accepted = FALSE)
  
  # If "synonym" is not null but no other sci names are given, "synonym" is actually
  # the scientific name of the species
  if(!is.null(synonyms) & is.null(sci_names)) 
    acc_sci_names_dat <- tibble(scientific_name = synonyms) %>%
    mutate(species = str_extract(scientific_name, accepted_species)) %>%
    filter(!is.na(species)) %>%
    mutate(accepted = TRUE)
  
  # Combine scientific names of synonyms and accepted names
  combined_sci_names_dat <- bind_rows(syn_sci_names_dat, acc_sci_names_dat)
  
  # Join to accepted species with taxon ID
  if(nrow(combined_sci_names_dat) > 0)
    species_dat <- full_join(species_dat, combined_sci_names_dat, by = c("species", "accepted")) %>%
    fill(taxid)
  
  species_dat
  
}

#' Fetch taxonomic data from the NCBI taxonomy database for a set of taxon ids 
#'
#' @param tax_ids Character vector of taxon IDs
#' @param chunk_size Number of chunks to split the tax_ids into when querying
#' the NCBI database. Should be 200 or less. 200 is probably fine.
#'
#' @return Dataframe with columns "taxid", "species", "scientific_name", and "accepted"
#' @examples
#' fetch_taxonomy(c("2726184", "2605333", "872590"))
fetch_taxonomy <- function(tax_ids, chunk_size = 200) {
  
  # Make sure taxonomic IDs don't include any missing values
  tax_ids <- as.character(tax_ids)
  assertthat::assert_that(is.character(tax_ids))
  assertthat::assert_that(!any(is.na(tax_ids)))
  
  if(length(tax_ids) < chunk_size) {
    fetch_taxonomy_chunk(tax_ids) %>%
      map_df(parse_tax_list)
  } else {
    n <- length(tax_ids)
    r <- rep(1:ceiling(n/chunk_size), each = chunk_size)[1:n]
    split(tax_ids, r) %>% 
      map(fetch_taxonomy_chunk) %>%
      # Remove one level of hierarchy to combine lists
      purrr::flatten() %>%
      map_df(parse_tax_list)
  }
  
}

#' Clean up species names downloaded from NCBI taxonomy database
#'
#' @param ncbi_names_raw Tibble with columns `taxid` `species` `accepted` 
#' and `scientific_name`. Names downloaded from NCBI taxonomy database
#' with fetch_taxonomy()
#'
#' @return Tibble; names with duplicates removed
#' 
clean_ncbi_names <- function(ncbi_names_raw) {
  # Clean up ncbi_names: some species have multiple accepted names.
  # These seem to include the species name w/o author and with author
  ncbi_accepted_mult_fixed <- ncbi_names_raw %>%
    filter(accepted == TRUE) %>%
    filter(!is.na(scientific_name)) %>%
    # Assume the sci. name. with most spaces has the author, and we want it
    add_count(taxid) %>%
    filter(n > 1) %>%
    select(-n) %>%
    mutate(n_spaces = str_count(scientific_name, " ")) %>%
    group_by(taxid) %>%
    # Discard ties (may lose some candidate names here, but so be it)
    slice_max(order_by = n_spaces, n = 1, with_ties = FALSE) %>%
    select(-n_spaces) 
  
  # Remove problematic names from original data, then add fixed names
  ncbi_names_raw %>%
    anti_join(ncbi_accepted_mult_fixed, by = "taxid") %>%
    bind_rows(ncbi_accepted_mult_fixed) %>%
    mutate(
      # Remove brackets around species name
      # (notation in NCBI taxonomic db that genus level taxonomy is uncertain)
      species = str_remove_all(species, "\\[|\\]"),
      # Remove year after authorship in scientific name
      scientific_name = str_remove_all(scientific_name, ", [0-9][0-9][0-9][0-9]")
    )
  
}

#' Exclude invalid names from taxonomic name resolution
#'
#' @param ncbi_names Names downloaded from NCBI taxonomy database.
#' Output of clean_ncbi_names()
#'
#' @return Dataframe
#' 
exclude_invalid_ncbi_names <- function(ncbi_names) {
  # Exclude names from consideration that aren't fully identified to species, or are hybrids
  # (assume there should be at least one space between genus and species)
  ncbi_names_exclude <-
    ncbi_names %>%
    filter(str_detect(species, " sp\\.| aff\\.| cf\\.| x ") | str_detect(scientific_name, " sp\\.| aff\\.| cf\\.| x ") | str_count(species, " ") < 1 | str_count(scientific_name, " ") < 1)
  
  ncbi_names %>%
    anti_join(ncbi_names_exclude, by = "species", na_matches = "never") %>%
    anti_join(ncbi_names_exclude, by = "scientific_name", na_matches = "never")
}

#' Select NCBI names for first round of taxonomic name resolution
#' 
#' Names that are considered "accepted" by NCBI and have full scientific name (with author)
#'
#' @param ncbi_names Names downloaded from NCBI taxonomy database.
#'
#' @return Dataframe (tibble)
#' 
select_ncbi_names_round_1 <- function(ncbi_names) {
  ncbi_names %>%
    filter(accepted == TRUE) %>%
    filter(!is.na(scientific_name)) %>%
    assertr::assert(is_uniq, taxid, scientific_name)
}


#' Classify results of taxon-tools matching
#'
#' @param match_results Dataframe; output of tt_match_names()
#'
#' @return Dataframe with column `result_type` added
#' 
tt_classify_result <- function(match_results) {
  match_results %>%
    add_count(query) %>%
    mutate(
      result_type = case_when(
        match_type != "no_match" & n == 1 ~ "single_match",
        match_type != "no_match" & n > 1 ~ "mult_match",
        match_type == "no_match" ~ "no_match",
        TRUE ~ NA_character_
      )
    ) %>%
    assert(not_na, result_type) %>%
    select(-n)
}

#' Resolve synonyms after matching names
#'
#' @param match_results_classified Dataframe; output of tt_classify_result() 
#' @param taxonomy_data Dataframe; taxonomic data with columns 
#' `taxon_id` `accepted_name_taxon_id` `status` `scientific_name`
#'
#' @return Dataframe; match results with column name `accepted_name` added.
#' Failures have `NA` for `accepted_name`.
#' 
tt_resolve_synonyms <- function(match_results_classified, taxonomy_data) {
  match_results_classified_with_taxonomy <-
    match_results_classified %>% 
    select(query, reference, match_type, result_type) %>%
    left_join(taxonomy_data, by = c(reference = "scientific_name"))
  
  accepted_single_match <-
    match_results_classified_with_taxonomy %>%
    filter(status == "accepted" & result_type == "single_match") %>%
    mutate(accepted_name = reference) %>%
    select(query, reference, accepted_name, match_type, status)
  
  accepted_single_synonyms <-
    match_results_classified_with_taxonomy %>% 
    filter(status == "synonym") %>%
    left_join(
      select(taxonomy_data, taxon_id, accepted_name = scientific_name), 
      by = c(accepted_name_taxon_id = "taxon_id"))  %>%
    select(query, reference, accepted_name, match_type, status) %>%
    group_by(query) %>%
    # Add count of number of resolved, accepted names per query
    mutate(n = n_distinct(accepted_name)) %>%
    ungroup() %>%
    filter(n == 1) %>%
    select(-n)
  
  success <- bind_rows(accepted_single_match, accepted_single_synonyms)
  
  failure <-
    match_results_classified_with_taxonomy %>%
    select(query, reference, match_type, status) %>%
    anti_join(success, by = "query")
  
  bind_rows(success, failure) %>% 
    verify(all(query %in% match_results_classified$query)) %>%
    verify(all(match_results_classified$query %in% query))
}

#' Select NCBI names for second round of taxonomic name resolution
#' 
#' Names that are considered synonyms by NCBI and have full scientific name (with author)
#'
#' @param match_results_resolved Dataframe; output of tt_resolve_synonyms() 
#' @param ncbi_names Names downloaded from NCBI taxonomy database.
#'
#' @return Dataframe (tibble)
#' 
select_ncbi_names_round_2 <- function(match_results_resolved_round_1, ncbi_names) {
  
  # Get IDs of all resolved names from round 1
  ncbi_id_resolved <-
    match_results_resolved_round_1 %>%
    left_join(ncbi_names, by = c(query = "scientific_name")) %>%
    filter(!is.na(accepted_name)) %>%
    assert(not_na, taxid)
  
  # Filter query names to those that failed round 1,
  # then to synonyms with a scientific name
  ncbi_names %>%
    anti_join(ncbi_id_resolved, by = "taxid") %>%
    filter(accepted == FALSE, !is.na(scientific_name)) %>%
    assertr::assert(is_uniq, scientific_name)
  
}

#' Select NCBI names for third round of taxonomic name resolution
#' 
#' NCBI names that are species only (lack scientific name)
#'
#' @param match_results_resolved_round_1 Dataframe; output of tt_resolve_synonyms() 
#' @param match_results_resolved_round_2 Dataframe; output of tt_resolve_synonyms() 
#' @param ncbi_names Names downloaded from NCBI taxonomy database.
#'
#' @return Dataframe (tibble)
#' 
select_ncbi_names_round_3 <- function(match_results_resolved_round_1, match_results_resolved_round_2, ncbi_names) {
  
  # Get IDs of all resolved names from round 1
  ncbi_id_resolved <-
    bind_rows(
      match_results_resolved_round_1,
      match_results_resolved_round_2) %>%
    left_join(ncbi_names, by = c(query = "scientific_name")) %>%
    filter(!is.na(accepted_name)) %>%
    assert(not_na, taxid)
  
  # Filter query names to those that failed round 1 + 2,
  # then to those lacking a scientific name (species name only)
  ncbi_names %>%
    anti_join(ncbi_id_resolved, by = "taxid") %>%
    filter(is.na(scientific_name)) %>%
    assertr::assert(is_uniq, species)
  
}

#' Combine results of name matching rounds 1-3
#'
#' @param ncbi_names_query Names downloaded from NCBI taxonomy database. 
#' @param ... Dataframes; output of tt_resolve_synonyms() 
#'
#' @return Dataframe; match results with NCBI taxid added
#'
combined_match_results <- function(ncbi_names_query, ...) {
  bind_rows(...) %>%
    left_join(select(ncbi_names_query, taxid, scientific_name), by = c(query = "scientific_name")) %>%
    left_join(select(ncbi_names_query, taxid, species), by = c(query = "species")) %>%
    mutate(taxid = coalesce(taxid.x, taxid.y)) %>%
    assert(not_na, taxid) %>%
    select(-taxid.x, -taxid.y)
}

#' Map resolved, accepted names to NCBI names
#'
#' @param match_results_resolved_all Dataframe; output of combined_match_results()
#'
#' @return Dataframe; NCBI names mapped to the accepted name in World Ferns
#' Does not include names that could not be matched to a single accepted name
#' @export
#'
#' @examples
make_ncbi_accepted_names_map <- function(match_results_resolved_all) {
  match_results_resolved_all %>%
    filter(!is.na(accepted_name)) %>% 
    select(taxid, accepted_name) %>%
    unique() %>% 
    assert(is_uniq, taxid) %>%
    # Add taxon (e.g., 'Foogenus barspecies fooinfraspname')
    mutate(
      rgnparser::gn_parse_tidy(accepted_name) %>% 
        select(taxon = canonicalsimple)
    ) %>%
    mutate(taxon = str_replace_all(taxon, " ", "_"))
}