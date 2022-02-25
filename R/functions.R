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
      taxonID = dwc_taxon_id, 
      acceptedNameUsageID = dwc_accepted_name_usage_id, 
      taxonomicStatus = dwc_taxonomic_status,
      rank = dwc_taxon_rank, 
      scientificName = dwc_scientific_name) %>%
    # Keep only species level and below
    filter(rank %in% c("form", "infraspecific name", "species", "subform", "subspecies", "subvariety", "variety")) %>%
    # Filter some names that were incorrectly labeled species level
    filter(str_detect(scientificName, "Polypodiaceae tribe Thelypterideae|Asplenium grex Triblemma|Pteridaceae tribus Platyzomateae|Filicaceae tribus Taenitideae", negate = TRUE)) %>%
    select(-rank) %>%
    # Note: taxonID is unique, but scientificName may not be (esp in case of ambiguous synonyms)
    assert(not_na, taxonID, scientificName) %>% 
    assert(is_uniq, taxonID)
  
  # Make sure all synonyms map correctly
  fow_accepted <- 
    fow %>%
    filter(str_detect(taxonomicStatus, "accepted")) %>%
    select(taxonID, scientificName, -taxonomicStatus)
  
  fow_synonyms <- 
    fow %>%
    filter(str_detect(taxonomicStatus, "synonym")) %>%
    select(taxonID, acceptedNameUsageID, scientificName, -taxonomicStatus)
  
  fow_synonyms %>%
    anti_join(fow_accepted, by = c(acceptedNameUsageID = "taxonID")) %>%
    verify(nrow(.) == 0, success_fun = success_logical)
  
  # Make sure all accepted names and synonyms are accounted for
  bind_rows(fow_accepted, fow_synonyms) %>%
    assert(is_uniq, taxonID) %>%
    anti_join(fow, by = "taxonID") %>%
    verify(nrow(.) == 0, success_fun = success_logical)
  
  # FIXME: Make sure accepted names and synonyms are distinct
  # A few repeats. Leave these in for now, but will need to fix.
  # fow_accepted %>%
  #   inner_join(fow_synonyms, by = c(name = "synonym")) %>%
  #   select(scientific_name = name) %>% 
  #   left_join(fow)
  
  fow
}

#' Load reference alignment files
#'
#' @param ref_aln_files Path to alignment (fasta) files,
#' named like 'ref_aln_atpA.fasta', where "atpA" is the target (gene) name
#' @return Tibble with two columns: "target" with name of gene, 
#' and "align_trimmed" with DNA sequences in a list
load_ref_aln <- function(ref_aln_files) {
  tibble(path = ref_aln_files) %>%
    mutate(
      target = str_match(path, "([a-zA-Z0-9\\-]+)\\.fasta") %>% 
        magrittr::extract(,2),
      align_trimmed = map(path, ~ape::read.FASTA(.) %>% as.matrix)
    ) %>%
    select(target, align_trimmed) %>%
    assert(not_na, everything()) %>%
    assert(is_uniq, target)
}

# Download Sanger sequences from GenBank----

#' Format a query string to download fern sequences from GenBank
#'
#' Helper function for fetch_fern_ref_seqs()
#'
#' @param target String; name of spacer region
#' @param start_date String; start date to include sequences
#' @param end_date String; end date to include sequences
#' @param strict Logical; should the search be strict with regards to gene name?
#'
#' @return String
format_fern_query <- function(target, start_date = "1980/01/01", end_date, strict = FALSE) {
  
  # Format query
	# Assume that we only want single genes or small sets of genes, not entire plastome.
	# Set upper limit to 7000 bp (we will fetch plastomes >7000 bp separately).

  # Define query based on whether it is a spacer region or not
  is_spacer <- str_detect(target, "-")
  if(isTRUE(is_spacer)) {
    # `spacer` should have exactly one hyphen
    assertthat::assert_that(isTRUE(str_count(target, "-") == 1))
    # split into names of two flanking genes
    flank_1 <- str_split(target, "-") %>%
      magrittr::extract2(1) %>%
      magrittr::extract2(1)
    flank_2 <- str_split(target, "-") %>%
      magrittr::extract2(1) %>%
      magrittr::extract2(2)
    if(isTRUE(strict)) {
      query <- glue('({flank_1}[GENE] AND {flank_2}[GENE]) AND Polypodiopsida[ORGN] AND 1:7000[SLEN] AND ("{start_date}"[PDAT]:"{end_date}"[PDAT]) NOT pseudogene') # nolint
    } else {
      query <- glue('({flank_1} OR {flank_2}) AND Polypodiopsida[ORGN] AND 1:7000[SLEN] AND ("{start_date}"[PDAT]:"{end_date}"[PDAT]) NOT pseudogene') # nolint
    }
  } else {
    if(isTRUE(strict)) {
      query <- glue('{target}[GENE] AND Polypodiopsida[ORGN] AND 1:7000[SLEN] AND ("{start_date}"[PDAT]:"{end_date}"[PDAT]) NOT pseudogene') # nolint
    } else {
      # some rbcL sequences use
      # "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit"
      # but not "rbcL"
      if(str_to_lower(target) == "rbcl") target <- "(ribu* OR rbcL)"
      query <- glue('{target} AND Polypodiopsida[ORGN] AND 1:7000[SLEN] AND ("{start_date}"[PDAT]:"{end_date}"[PDAT]) NOT pseudogene') # nolint
    }
  }
  query
}

#' Fetch raw sequences from GenBank
#'
#' Helper function for fetch_fern_ref_seqs()
#'
#' If ret_type is "gb", results will be in GenBank flat file format, as a
#' single string. Records are separated by '\\' within the string.
#'
#' @param query String; query string (formatted as if querying https://www.ncbi.nlm.nih.gov/nuccore/)
#' @param ret_type Return type; choose from "gb" 
#'   (raw text in GenBank flat-file format) or "fasta" (list of class DNAbin)
#' @param clean_names Logical; should names be cleaned to accession number only?
#'
#' @return String or list of class DNAbin
#' 
fetch_gb_raw <- function(query, ret_type = c("gb", "fasta"), clean_names = TRUE) {
  
  assertthat::assert_that(assertthat::is.string(query))
  assertthat::assert_that(assertthat::is.string(ret_type))
  assertthat::assert_that(assertthat::is.flag(clean_names))
  
  # Get list of GenBank IDs (GIs)
  uid <- reutils::esearch(term = query, db = "nucleotide", usehistory = TRUE)
  
  # Extract number of hits and print
  num_hits <- reutils::content(uid, as = "text") %>% str_match("<eSearchResult><Count>([:digit:]+)<\\/Count>") %>% magrittr::extract(,2)
  message(glue("Found {num_hits} accessions (UIDs) for query '{query}'"))
  
  # Download complete GenBank record for each and write it to a temporary file
  temp_dir <- tempdir()
  temp_file <- fs::path(temp_dir, digest::digest(query))
  
  # Make sure temp file doesn't already exist
  if (fs::file_exists(temp_file)) fs::file_delete(temp_file)
  
  # Download data
  reutils::efetch(uid, "nucleotide", rettype = ret_type, retmode = "text", outfile = temp_file)
  
  # Load results back in to R
  if (ret_type == "fasta") {
    results <- ape::read.FASTA(temp_file)
    # Clean names: keep only accession number
    if (isTRUE(clean_names)) {
      names(results) <- str_split(names(results), " ") %>% map_chr(1)
      names(results) <- str_remove_all(names(results), "\\.[0-9]$")
    }
  } else if (ret_type == "gb") {
    results <- readr::read_file(temp_file)
  } else {
    error("'ret_type' must be 'gb' or 'fasta'")
  }
  
  # Clean up
  fs::file_delete(temp_file)
  
  results
  
}

#' Parse a GenBank flatfile and retrieve integeneic spacer sequences
#'
#' Helper function for fetch_fern_ref_seqs()
#'
#' @param gb_raw String (character vector of length 1);
#' GenBank flatfile, possibly including multiple entries
#' @param target String; name of spacer to extract (e.g., "trnL-trnF")
#' @param req_intron Logical; should the intron in the left-side
#' flanking gene be included?
#' @param workers  Number of CPUs to run in parallel during parsing
#'
#' @return Tibble including columns for "seq", "error", and "target"
parse_gb_spacer <- function(gb_raw, target, req_intron = FALSE, workers) {

  # Change back to sequential when done (including on failure)
  on.exit(future::plan(future::sequential), add = TRUE)
  
  assertthat::assert_that(assertthat::is.string(gb_raw))
  assertthat::assert_that(assertthat::is.number(workers))
  
  # Set backend for parallelization
  future::plan(future::multisession, workers = workers)
  
  seqs_and_errors <-
    # Start with raw GenBank flat-file as single text string
    gb_raw %>%
    # '\\' is delimiter between entries; split up into one string each
    stringr::str_split("\n\\/\\/\n") %>%
    unlist %>%
    # Drop the last item, as it is just an empty line (after the last '\\')
    magrittr::extract(-length(.)) %>%
    # Can't start with any empty lines
    map(~str_remove(., "^[\n]+")) %>%
    # Loop over entries and extract spacer
    # genbankr::readGenBank() is slow, so do in parallel
    furrr::future_map(
      ~extract_spacer(
        ., target = target, req_intron = req_intron)) %>%
    transpose()

  # Close parallel workers
  future::plan(future::sequential)
  
  # Convert results to tibble
  seqs <- seqs_and_errors[["result"]] %>%
    purrr::compact() %>%
    do.call(c, .) %>%
    dnabin_to_seqtbl() %>%
    separate(accession, c("accession", "species"), sep = "__", fill = "right")

  tibble(
    seq = list(seqs), # seq is a nested tibble
    error = seqs_and_errors["error"], # error is a list-column with errors
    target = target
  )

}

#' Parse a GenBank flatfile and retrieve gene sequences
#'
#' Helper function for fetch_fern_ref_seqs()
#'
#' @param gb_raw String (character vector of length 1);
#' GenBank flatfile, possibly including multiple entries
#' @param target String; name of gene to extract
#' @param workers Number of CPUs to run in parallel during parsing
#'
#' @return Tibble including columns for "seq", "error", and "target"
parse_gb_gene <- function(gb_raw, target, workers) {

  # Change back to sequential when done (including on failure)
  on.exit(future::plan(future::sequential), add = TRUE)
  
  assertthat::assert_that(assertthat::is.string(gb_raw))
  assertthat::assert_that(assertthat::is.number(workers))
  
  # Set backend for parallelization
  future::plan(future::multisession, workers = workers)
  
  seqs_and_errors <-
    # Start with raw GenBank flat-file as single text string
    gb_raw %>%
    # '\\' is delimiter between entries; split up into one string each
    stringr::str_split("\n\\/\\/\n") %>%
    unlist %>%
    # Drop the last item, as it is just an empty line (after the last '\\')
    magrittr::extract(-length(.)) %>%
    # Can't start with any empty lines
    map(~str_remove(., "^[\n]+")) %>%
    # Loop over entries and extract gene
    # genbankr::readGenBank() is slow, so do in parallel
    furrr::future_map(~extract_gene(., target = target)) %>%
    transpose()

  # Close parallel workers
  future::plan(future::sequential)
  
  # Convert results to tibble
  seqs <- seqs_and_errors[["result"]] %>%
    purrr::compact() %>%
    do.call(c, .) %>%
    dnabin_to_seqtbl()
    
  # Split accession and species columns
  if(nrow(seqs > 0)) {
  seqs <-
    seqs %>%
    separate(accession, c("accession", "species"), sep = "__", fill = "right")
  }

  tibble(
    seq = list(seqs), # seq is a nested tibble
    error = seqs_and_errors["error"], # error is a list-column with errors
    target = target
  )

}

#' Extract a DNA sequence for a single spacer region from a genbank flat file
#'
#' Helper function for parse_gb_spacer()
#'
#' @param gb_entry String (character vector of length 1);
#' single entry from a genbank flatfile
#' @param target String; name of genes containing the spacer.
#' Must be formatted as two gene names separated by hyphen; e.g., 'trnL-trnF'
#' @param req_intron Logical; should an intron in the left-flanking gene
#' (in the case of trnL-trnF, trnL) be included?
#'
#' @return DNA sequence
#' 
extract_spacer_fragile <- function (gb_entry, target, req_intron = FALSE) {

  # `target` should have exactly one hyphen
	assertthat::assert_that(isTRUE(str_count(target, "-") == 1))
	
	# split into names of two flanking genes
	flank_1 <- str_split(target, "-") %>% magrittr::extract2(1) %>% magrittr::extract2(1)
	flank_2 <- str_split(target, "-") %>% magrittr::extract2(1) %>% magrittr::extract2(2)

  # Parse the GenBank file
  gb_parsed <- suppressMessages(genbankr::readGenBank(text = gb_entry, partial = TRUE, verbose = FALSE))

  if(is.null(gb_parsed)) {
    stop("Could not parse GenBank file")
  }
  
  # Extract the accession
  accession <- attributes(gb_parsed)$accession
  if(is.null(accession)) {
    stop(glue::glue("No accession detected"))
  }
  # Extract the feature locations
  other_features <-
    genbankr::otherFeatures(gb_parsed) %>% 
    as_tibble() %>%
    mutate(across(where(is.factor), as.character))

  if(!"type" %in% colnames(other_features)) {
    stop(glue::glue("No 'type' for {target} gene detected in accession {accession}"))
  }

  # Optionally check for intron in flank_1  
  if(isTRUE(req_intron)) {
    # Return NULL if no gene detected
    if(!"gene" %in% colnames(other_features)) {
      stop(glue::glue("No {flank_1} gene detected in accession {accession}"))
    }
    # Find position of intron
    intron_features <-
      other_features %>%
      filter(type == "intron", str_detect(gene, regex(flank_1, ignore_case = TRUE)))
    # Error if no intron detected
    if(nrow(intron_features) == 0) {
      stop(glue::glue("No intron for {flank_1} gene detected in accession {accession}"))
    }
    if(nrow(intron_features) > 1) {
      stop(glue::glue("Multiple introns for {flank_1} gene detected in accession {accession}"))
    }
  }
  # Error if no note detected
  if(!"note" %in% colnames(other_features)) {
    stop(glue::glue("No 'note' detected in accession {accession}"))
  }

  # Find position of spacer
  spacer_features <-
    other_features %>%
    filter(
      type == "misc_feature", 
      str_detect(note, regex(as.character(glue::glue("{flank_1}-{flank_2}")), ignore_case = TRUE)),
      str_detect(note, "spacer"),
      str_detect(note, "intron", negate = TRUE)
      )
   # Error if no spacer detected
  if(nrow(spacer_features) == 0) {
    stop(glue::glue("No spacer for {target} gene detected in accession {accession}"))
  }
  if(nrow(spacer_features) > 1) {
    stop(glue::glue("Multiple spacer for {target} gene detected in accession {accession}"))
  }
  # Error if start/end don't agree
  if(isTRUE(req_intron)) {
    if(intron_features$start >= spacer_features$end) {
      stop(glue::glue("Start intron after spacer end for {target} gene detected in accession {accession}"))
    }
  }

  # If no errors, extract sequence
  if(isTRUE(req_intron)) {
    # Extract the intron + spacer sequence, or
    seq <-
      Biostrings::getSeq(gb_parsed) %>% 
      Biostrings::subseq(intron_features$start, spacer_features$end) %>%
      ape::as.DNAbin()
    } else {
    # extract the spacer sequence
    seq <-
      Biostrings::getSeq(gb_parsed) %>% 
      Biostrings::subseq(spacer_features$start, spacer_features$end) %>%
      ape::as.DNAbin()
  }

  # Name as the sequence as accession__organism
  organism <- attributes(gb_parsed)$sources$organism

  # attributes(gb_parsed)$sources$organism
  names(seq) <- glue::glue("{accession}__{organism}") %>% as.character()
  
  seq
  
}

extract_spacer <- purrr::safely(extract_spacer_fragile)

#' Extract a DNA sequence for a single gene from a genbank flat file
#'
#' Helper function for parse_gb_gene()
#'
#' @param gb_entry String (character vector of length 1);
#' single entry from a genbank flatfile
#' @param target String; name of gene
#'
#' @return DNA sequence
#' 
extract_gene_fragile <- function (gb_entry, target) {

  # Parse the GenBank file
  gb_parsed <- suppressMessages(genbankr::readGenBank(text = gb_entry, partial = TRUE, verbose = FALSE))

  if(is.null(gb_parsed)) {
    stop("Could not parse GenBank file")
  }
  
  # Extract the accession
  accession <- attributes(gb_parsed)$accession
  if(is.null(accession)) {
    stop(glue::glue("No accession detected"))
  }

  # Extract genes
  genes <- genbankr::genes(gb_parsed) %>% as_tibble()

  assertthat::assert_that(
    nrow(genes) > 0,
    msg = glue("No genes detected in accession {accession}"))

  target_genes <- filter(genes, gene == target)

  assertthat::assert_that(
    nrow(target_genes) > 0,
    msg = glue("Target gene {target} not detected in accession {accession}"))

  assertthat::assert_that(
    is.numeric(target_genes$start))
  
  assertthat::assert_that(
    is.numeric(target_genes$end))
  
  assertthat::assert_that(target_genes$start < target_genes$end)

  # Extract the gene
  seq <-
    Biostrings::getSeq(gb_parsed) %>% 
    Biostrings::subseq(target_genes$start, target_genes$end) %>%
    ape::as.DNAbin()

  # Name as the sequence as accession__organism
  organism <- attributes(gb_parsed)$sources$organism

  # attributes(gb_parsed)$sources$organism
  names(seq) <- glue::glue("{accession}__{organism}") %>% as.character()
  
  seq
  
}

extract_gene <- safely(extract_gene_fragile)

#' Download fern DNA sequences to use as reference for extracting
#' target sequences with BLAST
#'
#' @param target Name of target locus (e.g, "rbcL", "trnL-trnF")
#' @param start_date Start date for filtering GenBank accessions
#' @param end_date End date for filtering GenBank accessions
#' @param req_intron Logical; should the intronic region be included?
#' Only applies to trnL-trnF
#' @param workers Number of CPUs to run in parallel during parsing of GenBank flatfile
#' @param strict Logical; use strict version of search string or not?
#' @param accs_exclude Character vector of GenBank accessions to exclude
#' from results
#'
#' @return Tibble of DNA sequences
fetch_fern_ref_seqs <- function(target, start_date = "1980/01/01", end_date, 
  req_intron = FALSE, workers, strict = FALSE, accs_exclude = NULL) {

   # Check if query target is a spacer region
  is_spacer <- str_detect(target, "-")

  # Download sequences in GenBank flatfile format (plain text)
  gb_raw <-
  # Format query
  format_fern_query(
    target = target,
    start_date = start_date,
    end_date = end_date,
    strict = strict) %>%
  # Fetch raw GenBank flatfile
  fetch_gb_raw(query = ., ret_type = "gb")

  # Parse features
  # Result will be sequences in a tibble
  if(is_spacer) {
    res <- parse_gb_spacer(
        gb_raw = gb_raw, target = target, req_intron = req_intron,
        workers = workers
    )
  } else {
    res <- parse_gb_gene(
        gb_raw = gb_raw, target = target,
        workers = workers
    )
  }

  # Remove accessions in exclude list
  if (!is.null(accs_exclude)) {
    res$seq[[1]] <- filter(res$seq[[1]], !accession %in% accs_exclude)
  }

  res

}

#' Filter a set of sequences
#'
#' Selects the sequence with most non-missing bases per genus per target locus
#'
#' @param fern_ref_seqs_raw Output of fetch_fern_ref_seqs()
#' @return Tibble of filtered sequences
filter_ref_seqs <- function(fern_ref_seqs_raw) {
  
# Filter based on non-missing bases
  fern_ref_seqs_raw %>%
    select(seq, target) %>%
    unnest(cols = c(seq, target)) %>%
    mutate(
      # Here, seq_len only includes a,c,t,g (no missing or ambiguous bases)
      seq_len = map_dbl(seq, count_actg),
      # remove brackets fromm species names
      species = str_remove_all(species, "\\[|\\]"),
      # extract genus
      genus = str_split(species, " ") %>% map_chr(1)) %>%
    # filter by most non-missing bases per genus per target
    group_by(genus, target) %>%
    slice_max(order_by = seq_len, with_ties = FALSE) %>%
    ungroup()
}

#' Download a set of fern sequences for a given target region
#'
#' @param target Name of target
#' @param start_date Earliest date to download
#' @param end_date Most recent date to download
#' @param accs_exclude Character vector of accessions to exclude
#' @param terms_remove GREP string to remove from `target` when querying
#' @param strict Logical; should the search be strict with regards to gene name?
#'
#' @return Dataframe with one row per sequence and columns
#' for the sequence and accession
#'
#' fetch_fern_sanger_seqs("rbcL", start_date = "2018/01/01", end_date = "2018/01/10")
fetch_fern_sanger_seqs <- function(
  target, start_date = "1980/01/01", end_date, accs_exclude = NULL,
  strict = FALSE) {

  assertthat::assert_that(assertthat::is.string(target))
  assertthat::assert_that(assertthat::is.string(end_date))
  assertthat::assert_that(assertthat::is.string(start_date))

 seqs <-
  # Format query
  format_fern_query(
    target = target,
    start_date = start_date,
    end_date = end_date,
    strict = strict) %>%
  # Download sequences
  fetch_gb_raw(ret_type = "fasta")

  # Exclude any accessions in exclusion list
  if (!is.null(accs_exclude)) seqs <- seqs[!seqs %in% accs_exclude]

  # Convert to tibble
  tibble::tibble(
    seq = split(seqs, seq_along(seqs)), accession = names(seqs), gene = target)

}

#' Extract sequences by querying against a reference library
#'
#' Wrapper for supercrunch Reference_Blast_Extract.py module
#'
#' @param query_seqtbl DNA sequences to extract, formatted as sqtbl
#' @param ref_seqtbl DNA sequences to use as a BLAST reference database, formatted as
#' tibble with list-column containing reference alignments called "align_trimmed"
#' @param target Name of target locus
#' @param blast_flavor Name of BLAST algorithm to use
#' @param other_args Additional arguments formatted as a character vector to send
#' to superCRUNCH Reference_Blast_Extract.py
#' @param echo Logical; should the output of Reference_Blast_Extract.py be printed to the screen?
#' @param blast_res Logical; should the BLAST output be included in the results?
#'
#' @return Tibble

extract_from_ref_blast <- function(query_seqtbl, ref_seqtbl, target, blast_flavor, other_args = NULL, echo = FALSE, blast_res = FALSE) {
  
  # To avoid confusion with columns named 'target'
  target_select <- target

  # Filter sequences, convert to DNAbin:
  # - query (filter on target only)
  query_seqs <- query_seqtbl %>%
    filter(gene == target_select) %>%
    unique() %>%
    verify(nrow(.) > 0, error_fun = err_msg("No query sequences matching target")) %>%
    seqtbl_to_dnabin(name_col = "accession", seq_col = "seq") %>%
    # remove gaps
    ape::del.gaps()
  
  # - reference (filter on target, and optionally by feature)
  ref_seqs <- ref_seqtbl %>%
    filter(target == target_select) %>%
    verify(nrow(.) > 0, error_fun = err_msg("No reference sequences matching target")) %>%
    # Assuming sequences are in list-col "align_trimmed"
    pull(align_trimmed) %>%
    magrittr::extract2(1) %>%
    # remove gaps (important for ref database)
    ape::del.gaps()
  
  # Create temporary input and output folders
  sha <- digest::digest(c(query_seqs, ref_seqs))
  in_folder <- tempfile(glue("{sha}_in"))
  out_folder <- tempfile(glue("{sha}_out"))
  
  if(fs::dir_exists(in_folder)) fs::dir_delete(in_folder)
  if(fs::dir_exists(out_folder)) fs::dir_delete(out_folder)
  
  fs::dir_create(in_folder)
  fs::dir_create(out_folder)
  
  # Write out sequences
  query_file <- fs::path_ext_set(target, "fasta")
  ref_file <- paste0(target, "_ref") %>% fs::path_ext_set("fasta")
  
  ape::write.FASTA(query_seqs, fs::path(in_folder, query_file))
  ape::write.FASTA(ref_seqs, fs::path(in_folder, ref_file))
  
  # Set up super-crunch arguments
  args <- c(
    "Reference_Blast_Extract.py",
    "-i", fs::path_abs(in_folder),
    "-o", fs::path_abs(out_folder),
    "-e", query_file,
    "-d", ref_file,
    "-b", blast_flavor,
    other_args
  )
  
  # Run supercrunch
  sc_res <- processx::run("supercrunch", args, echo = echo)
  
  # Read in output files, return as list
  out_files <- list.files(out_folder, full.names = TRUE, recursive = TRUE)
  
  blast_file <- out_files[str_detect(out_files, "blast_results.txt")]
  log_file <- out_files[str_detect(out_files, "Log_File_")]
  bad_seqs_file <- out_files[str_detect(out_files, "Log_BadSeqs_")]
  extracted_seqs_file <- out_files[str_detect(out_files, "_extracted.fasta")]
  
  blast_res <- NULL
  log_res <- NULL
  bad_seqs_res <- NULL
  extracted_seqs_res <- NULL
  
  fmt6_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  
  
  if(length(blast_file) > 0 && isTRUE(blast_res)) blast_res <- suppressMessages(
    {read_tsv(blast_file, col_names = fmt6_cols, col_types = "ccdddddddddd")})
  if(length(log_file) > 0) log_res <- suppressMessages(read_tsv(log_file))
  if(length(bad_seqs_file) > 0) bad_seqs_res <- ape::read.FASTA(bad_seqs_file) %>% dnabin_to_seqtbl()
  if(length(extracted_seqs_file) > 0) extracted_seqs_res <- ape::read.FASTA(extracted_seqs_file) %>% dnabin_to_seqtbl()
  
  fs::dir_delete(in_folder)
  fs::dir_delete(out_folder)
  
  tibble::tibble(
    blast = list(blast_res),
    log = list(log_res),
    bad_seqs = list(bad_seqs_res),
    extracted_seqs = list(extracted_seqs_res),
    target = target_select,
    blast_flavor = blast_flavor,
    stdout = list(sc_res$stdout),
    stderr = list(sc_res$stderr)
  )
  
}

# Helper function to make tibble including custom error message,
# accession, and gene. Uses `accession` and `gene` in parent environment.
error_tbl <- function(accession, gene, msg) {
  tibble(accession = accession, gene = gene, msg = msg)
}

#' Extract a DNA sequence from a single entry
#' in a genbank flatfile
#'
#' @param gb_entry String (character vector of length 1);
#' single entry from a genbank flatfile
#' @param gene Name of gene
#'
#' @return List with two items:
#'   - seq: Named character vector; DNA sequence
#'   - error: Tibble with error message, gene name, and accession
#' 
extract_sequence <- function (gb_entry, gene, accs_exclude = NULL) {

  # Check for accession number
  accession_detected <- str_detect(gb_entry, "ACCESSION")
	accession_detected_msg <- assertthat::validate_that(
		accession_detected,
		msg = glue::glue("Genbank flatfile not valid (missing ACCESSION); no sequence extracted")
	)
  if (!accession_detected) {
		message(accession_detected_msg)
		return(
      list(
        seq = NULL,
        error = tibble(gene = gene, msg = accession_detected_msg)
      )
    )
	}

   # Extract accession number
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

  # If exclusion list is present and it's on the list, skip it
  if(!is.null(accs_exclude) && accession %in% accs_exclude) return (list(seq = NULL, error = NULL))
	
	# Check for FEATURES and ORIGIN fields
  features_detected <- str_detect(gb_entry, "FEATURES")
	features_detected_msg <- assertthat::validate_that(
		features_detected,
		msg = glue::glue("Genbank flatfile not valid (missing FEATURES); no sequence extracted")
	)
  if(!features_detected) return(
    list(seq = NULL, error = error_tbl(accession = accession, gene = gene, msg = features_detected_msg))
  )

  origin_detected <- str_detect(gb_entry, "ORIGIN")
	origin_detected_msg <- assertthat::validate_that(
		origin_detected,
		msg = glue::glue("Genbank flatfile not valid (missing ORIGIN); no sequence extracted")
	)
  if(!origin_detected) return(
    list(seq = NULL, error = error_tbl(accession = accession, gene = gene, msg = origin_detected_msg))
  )
  
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
    unlist %>%
    str_squish()
  
  # Make sure target gene is detected
  gene_detected <- any(str_detect(gene_range_list, regex(gene, ignore_case = TRUE)))
  gene_detected_msg <- assertthat::validate_that(
    gene_detected,
    msg = glue::glue("Gene {gene} not detected in accession {accession}")
  )
  if(!gene_detected) return(
    list(seq = NULL, error = error_tbl(accession = accession, gene = gene, msg = gene_detected_msg))
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
  
  # Check for duplicated genes
  gene_not_duplicated <- length(unique(gene_range)) <= 2
	gene_not_duplicated_msg <- assertthat::validate_that(
		gene_not_duplicated,
		msg = glue::glue("Duplicate copies of {gene} gene detected in accession {accession}")
	)
  if(!gene_not_duplicated) return(
    list(seq = NULL, error = error_tbl(accession = accession, gene = gene, msg = gene_not_duplicated_msg))
  )
	
  # Check that full range of gene was detected
  gene_full <- length(gene_range) > 1 && 
    is.numeric(gene_range) && 
    !anyNA(gene_range) && 
    gene_range[1] <= gene_range[2] && 
    gene_range[2] >= gene_range[1]
  
 	gene_full_msg <- assertthat::validate_that(
		gene_full,
		msg = glue::glue("Full range of {gene} gene not detected in accession {accession}")
	)
	if(!gene_full) return(
    list(seq = NULL, error = error_tbl(accession = accession, gene = gene, msg = gene_full_msg))
  )
  
  # Extract sequence, subset to target gene
	sequence <- 
		gb_entry %>%
		paste(sep = "") %>%
		str_remove_all("\n") %>%
		str_remove_all('\"') %>%
    strex::str_after_last("ORIGIN") %>% 
		str_remove_all(" ") %>%
		str_remove_all("[0-9]") %>%
		substr(gene_range[1], gene_range[2])

  # Check for valid DNA sequences
  # - make a grep query that will hit any non-IUPAC character (in upper case)
  non_iupac <- "[^A^C^G^T^U^R^Y^S^W^K^M^B^D^H^V^N^\\.^\\-^\\?]"

  iupac_only <- str_detect(str_to_upper(sequence), non_iupac, negate = TRUE)

  iupac_only_msg <- assertthat::validate_that(
		iupac_only,
		msg = glue::glue("Non-IUPAC characters detected in {flank_1}-{flank_2} spacer of accession {accession}")
	)

	if(!iupac_only) return(
    list(seq = NULL, error = error_tbl(accession = accession, gene = gene, msg = iupac_only_msg))
  )
  
  # If pass all checks, return sequence with no error
  list(
    # convert sequence to ape DNAseq
    seq = set_names(sequence, accession),
    error = NULL
  ) 
}

#' Clean up results from extract_from_ref_blast()
#'
#' @param extract_from_ref_blast_res results from extract_from_ref_blast()
#' @param blast_flavor_select Type of blast algorithm to use
#' @results Tibble (seqtbl)
#'
clean_extract_res <- function(extract_from_ref_blast_res, blast_flavor_select) {
  extract_from_ref_blast_res %>%
  filter(blast_flavor == blast_flavor_select) %>%
  select(extracted_seqs, target) %>%
  unnest(extracted_seqs)
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

#' Download a set of fern sequence metadata for a given target
#'
#' @param target Name of target
#' @param start_date Earliest date to download
#' @param end_date Most recent date to download
#' @param accs_exclude Character vector of accessions to exclude
#' @param strict Logical; use strict version of search string or not?
#'
#' @return Datarame
#' fetch_fern_metadata("rbcL", start_date = "2018/01/01", end_date = "2018/02/01")
#' fetch_fern_metadata("trnL-trnF_intron", start_date = "2017/01/01", end_date = "2018/02/01")
fetch_fern_metadata <- function(
  target, start_date = "1980/01/01", end_date, accs_exclude = NULL, strict = FALSE) {
  
  assertthat::assert_that(assertthat::is.string(target))
  assertthat::assert_that(assertthat::is.string(end_date))
  assertthat::assert_that(assertthat::is.string(start_date))
  
  # Format query
  query <- format_fern_query(
    target = target,
    start_date = start_date, 
    end_date = end_date,
    strict = strict)

  # Fetch standard metadata
  metadata <- fetch_metadata(
    query = query,
    # don't fetch `slen` (length of accession; will calculate length of actual sequence later)
    col_select = c("gi", "caption", "taxid", "title", "subtype", "subname")) %>%
    rename(accession = caption) %>%
    # GenBank accession should be non-missing, unique
    assert(not_na, accession) %>%
    assert(is_uniq, accession)

  # Exclude accessions in accs_exclude
  if (!is.null(accs_exclude)) {
    metadata <-
      metadata %>%
      filter(!accession %in% accs_exclude)
  }

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
    mutate(target = target)
  
}

#' Align sequences used as reference
#' @param fern_ref_seqs Reference DNA sequences formatted as seqtbl
#' @param target_select String; name of target locus
#' @param n_threads Number of threads to use when aligning sequences
#' @return Tibble (seqtbl)
align_ref_seqs <- function(fern_ref_seqs, target_select, n_threads = 1) {
  fern_ref_seqs %>%
    assert(not_na, species, accession, target) %>%
    filter(target == target_select) %>%
    mutate(species = str_replace_all(species, " ", "_")) %>%
    # name sequences as [accession]__[species]
    unite("seq_name", c(accession, species, target), sep = "__") %>%
    seqtbl_to_dnabin(name_col = "seq_name") %>%
    ips::mafft(
      x = .,
      options = "--adjustdirection",
      exec = "/usr/bin/mafft",
      thread = n_threads) %>%
    remove_mafft_r() %>%
    dnabin_to_seqtbl()
}

#' Filter raw fasta sequences to one per genus
#'
#' Used in prep_ref_seqs plan
#'
#' @param raw_fasta Raw fasta sequences in tibble
#' @param raw_meta Raw metadata in tibble
#'
#' @return Tibble with one sequence (longest) per genus per target locus
filter_raw_fasta_by_genus <- function(raw_fasta, raw_meta) {

  # Obtain NCBI taxon names for accessions
  ncbi_tax_names <-
  raw_meta %>%
    pull(taxid) %>%
    unique() %>%
    taxize::ncbi_get_taxon_summary(id = .) %>%
    as_tibble()

  # Further format metadata with taxon names
  # (some infrasp, but treat as species)
  meta_with_ncbi_names <-
  raw_meta %>%
    inner_join(ncbi_tax_names, by = c(taxid = "uid")) %>%
    # only keep taxa at species level or below
    filter(rank %in% c("forma", "species", "subspecies", "varietas")) %>%
    # drop square brackets in names
    mutate(name = str_remove_all(name, "\\[|\\]")) %>%
    # split out genus
    mutate(
      genus = str_split(name, " ") %>% map_chr(1),
      sp_epithet = str_split(name, " ") %>% map_chr(2)
      ) %>%
    unite("species", c(genus, sp_epithet), na.rm = TRUE, remove = FALSE) %>%
    select(accession, genus, species)

  # filter by longest sequence per genus per target locus
  raw_fasta %>%
    inner_join(meta_with_ncbi_names, by = "accession") %>%
    mutate(seqln = map_dbl(seq, ~length(.[[1]]))) %>%
    group_by(target, genus) %>%
    slice_max(order_by = "seqln", n = 1, with_ties = FALSE) %>%
    ungroup()
}


#' Trim a DNA alignment by removing sequences prior to a motif
#' 
#' e.g., trim the portion of rps4-trnS prior to the stop codon of rps4,
#' so that only the rps4-trnS spacer is retained (not rps4 gene)
#' 
#' Used in prep_ref_seqs plan
#'
#' @param aln_tbl Tibble including DNA alignment as a list-column
#' @param aln_col Name of list-column with DNA alignment
#' @param target_select Name of target locus (should be match one of the
#' values of the "target" column in aln_tbl)
#' @param motif DNA motif to match (won't match on indels)
#' @param expect_start Approximate expected start range
#' @param expect_end Approximate expected end range
#'
#' @return Tibble including DNA alignment as a list-column. The alignment for the
#' selected gene will be trimmed to remove all bases prior to and including the motif
#' 
trim_align_by_motif <- function(
	aln_tbl, aln_col = "align_trimmed",
	target_select = "rps4-trnS", motif = "CTTAATGA",
	expect_start = 500,
	expect_end = 700) {
	
	# Extract alignment, preserving gaps
	aln <-
		aln_tbl %>%
		filter(target == target_select) %>%
		pull(.data[[aln_col]]) %>%
		magrittr::extract2(1) %>%
		DNAbin_to_DNAstringset(remove_gaps = FALSE)
	
	# Search for pattern matching end of region
	# vmatchPattern can't account for indels in DNAStringSet
	matches <- Biostrings::vmatchPattern(
    Biostrings::DNAString(motif), aln, with.indels = FALSE)
	
  # Convert vmatchPattern output to tibble
	group_names_tbl <- tibble(
		group_name = names(matches),
		group = 1:length(matches))
	
	hits_tbl <-
		as.data.frame(matches) %>%
		as_tibble() %>%
		select(-group_name) %>%
		left_join(group_names_tbl, by = "group") %>%
		# Run checks: should only be a single unique match location,
		# within expected range
		verify(n_distinct(.$end) == 1) %>%
		verify(all(.$start > expect_start)) %>%
		verify(all(.$end < expect_end))
	
	# Trim alignment
	# removes all sequences before (and including) the matched motif
	aln_clean <-
		Biostrings::DNAStringSet(aln, start = (unique(hits_tbl$end) + 1)) %>%
		as.DNAbin() %>%
		as.matrix()
	
	# Convert back to tbl
	aln_tbl %>%
		filter(target != target_select) %>%
		bind_rows(
			tibble(
				target = target_select,
				!!aln_col := list(aln_clean)
			)
		)
}

#' Write out DNA alignments from a tibble
#'
#' Used in prep_ref_seqs plan
#'
#' @param tbl Tibble with aligned DNA sequences. Will only write out first row.
#' @param aln_col Name of list-column with aligned DNA sequences
#' @param gene_col Name of column with gene (target locus) names 
#' @param dir Directory to write output
#' @param prefix Characters to add to file name preceding the target locus.
#' @param postfix Characters to add to file name following the target locus.
#'
#' @return Path to output
#'
write_fasta_from_tbl <- function(
  tbl, aln_col = "align_trimmed",
  gene_col = "target",
  dir = "intermediates/ref_seqs",
  prefix = "",
  postfix = "_ref_aln.fasta") {
  write_fasta_tar(
    tbl[[aln_col]][[1]],
    fs::path(dir, paste0(prefix, tbl[[gene_col]][[1]], postfix))
  )
}

#' Write out phylogenetic trees from a tibble
#'
#' Used in prep_ref_seqs plan
#'
#' @param tbl Tibble with phylogenetic tree. Will only write out first row.
#' @param tree_col Name of list-column with phylogenetic tree
#' @param gene_col Name of column with gene (target locus) names 
#' @param dir Directory to write output
#' @param prefix Characters to add to file name preceding the target locus.
#' @param postfix Characters to add to file name following the target locus.
#'
#' @return Path to output
#'
write_tree_from_tbl <- function(
  tbl, tree_col = "tree",
  gene_col = "target",
  dir = "intermediates/ref_seqs", 
  prefix = "",
  postfix = "_ref_phy.tree") {
  write_tree_tar(
    tbl[[tree_col]][[1]],
    fs::path(dir, paste0(prefix, tbl[[gene_col]][[1]], postfix))
  )
}

# Taxonomic name resolution ----

# Specify varieties to exclude from collapsing during taxonomic name resolution
define_varieties_to_keep = function() {
  c(
	  "Arachniodes simplicior var. simplicior",
	  "Asplenium attenuatum var. attenuatum",
	  "Asplenium bulbiferum subsp. bulbiferum",
	  "Asplenium obovatum subsp. obovatum",
	  "Asplenium trichomanes subsp. trichomanes",
	  "Botrychium lanceolatum subsp. lanceolatum",
	  "Equisetum ramosissimum subsp. ramosissimum",
	  "Pellaea mucronata subsp. mucronata",
	  "Pleopeltis polypodioides var. polypodioides",
	  "Polystichum polyblepharum var. polyblepharum",
    "Dicksonia lanata subsp. lanata",
    "Dryopteris simasakii var. simasakii",
    "Elaphoglossum peltatum f. peltatum",
    "Hypolepis rugosula subsp. rugosula"
  )
}

#' Inspect results of name matching with taxastand
#' 
#' Subsets to fuzzily matched names and names without a match, prepares for updating database
#'
#' @param match_results_resolved_all Dataframe (tibble); results of taxonomic name matching
#' with {taxastand}; output of combined_match_results()
#' @return Dataframe (tibble) of fuzzily matched results, with several columns
#' pre-populated with default values to be inspected and used for updating the
#' taxonomic reference database. Includes:
#' - query: Query string
#' - matched_name: Matched name
#' - matched_status: Status of matched name
#' - query_match_taxon_agree: Logical; do the canonical names
#'   match between query and matched_name?
#' - use_query_as_synonym: 1 to use query as synonym of matched name
#'   in next build of taxonomic reference
#' - use_query_as_accepted: 1 to use query as accepted name, with matched
#'   name as a synonym, in next build of taxonomic reference
#' - use_query_as_new: 1 to use query as new name in next build of
#'   taxonomic reference
#' - taxonomicStatus, namePublishedIn, nameAccordingTo, taxonRemarks:
#'   Darwin Core columns to use in next build of taxonomic reference
inspect_ts_results <- function(match_results_resolved_all) {
  # - Fuzzily-matched names
  match_results_resolved_all %>%
    filter(match_type == "auto_fuzzy") %>%
    select(query, matched_name, matched_status) %>%
    mutate(
      rgnparser::gn_parse_tidy(query) %>%
        select(query_taxon = canonicalsimple),
      rgnparser::gn_parse_tidy(matched_name) %>%
        select(matched_taxon = canonicalsimple)
    ) %>%
    mutate(
      query_taxon = str_replace_all(
        query_taxon, "Vandenboschia radicans type", "Vandenboschia radicans"),
      query_match_taxon_agree = query_taxon == matched_taxon,
      # The below values are defaults. Each row needs to be checked manually.
      use_query_as_synonym = case_when(
        query_match_taxon_agree == TRUE ~ 1,
        TRUE ~ 0),
      use_query_as_accepted = 0,
      use_query_as_new = 0,
      namePublishedIn	= NA,
      nameAccordingTo = NA,
      taxonRemarks = case_when(
        query_match_taxon_agree == TRUE ~ "author variant",
        query_match_taxon_agree == FALSE ~ "spelling mistake in species name"
      ),
      notes = NA
      ) %>%
    select(-query_taxon, -matched_taxon) %>%
    arrange(query_match_taxon_agree, query) %>%
    unique() %>%
    # - Non-matching names
    bind_rows(
      filter(match_results_resolved_all, match_type == "no_match") %>%
        select(query) %>%
        unique() %>%
        arrange(query)
    )
}

# Combine Sanger sequences data ----

# Count the number of non-missing bases in a DNA sequence
#' @param seq List of class DNAbin of length one.
count_non_missing <- function (seq) {
  bases <- ape::base.freq(seq, all = TRUE, freq = TRUE)
  sum(bases[!names(bases) %in% c("n", "N", "-", "?")])
}

#' Combine sanger sequence metadata with sequences, join to resolved names
#' and filter by sequence length and if name was resolved or not
#' 
#' Drops sequences with scientific names that could not be resolved,
#' adds `otu` column ({species}|{accession}|{gene})
#'
#' @param raw_meta Sanger sequence metadata; output of fetch_fern_metadata()
#' @param raw_fasta Sanger sequences; output of fetch_fern_gene()
#' @param ncbi_accepted_names_map Dataframe mapping NCBI taxid to accepted
#' species name; output of make_ncbi_accepted_names_map()
#' @param min_gene_len Number: minimum length (bp) required for Sanger genes
#' @param min_spacer_len Number: minimum length (bp) required for Sanger intergenic spacer regions
#'
#' @return Tibble with Sanger sequence metadata, sequences, and accepted name
#' 
combine_and_filter_sanger <- function(
  raw_meta, raw_fasta, ncbi_accepted_names_map,
  min_gene_len, min_spacer_len) {
  
  # Check that input names match arguments
  check_args(match.call())
  
  # Join metadata and fasta sequences
  raw_meta %>%
    left_join(raw_fasta, by = c("accession", "target")) %>%
    # Filter out null sequences
    filter(!map_lgl(seq, is.null)) %>%
    # Inner join to name resolution results: will drop un-resolved names
    inner_join(ncbi_accepted_names_map, by = "taxid") %>%
    # Calculate actual seq length
    # only count non-missing bases
    mutate(seq_len = map_dbl(seq, ~count_non_missing(.[1]))) %>%
    # Categorize gene type
    mutate(
      target_type = case_when(
        str_detect(target, "-") ~ "spacer",
        TRUE ~ "gene"
      )
    ) %>%
    assert(not_na, target_type) %>%
    # Filter by minimum seq. length
    filter((seq_len > min_gene_len & target_type == "gene") | (seq_len > min_spacer_len & target_type == "spacer")) %>%
    assert(not_na, accession, seq, resolved_name, species) %>%
    # Create OTU column for naming sequences as species|accession|target
    # - first make sure there are no spaces or `|` in species, accession, or target
    verify(all(str_detect(species, " ", negate = TRUE))) %>%
    verify(all(str_detect(accession, " ", negate = TRUE))) %>%
    verify(all(str_detect(target, " ", negate = TRUE))) %>%
    verify(all(str_detect(species, "\\|", negate = TRUE))) %>%
    verify(all(str_detect(accession, "\\|", negate = TRUE))) %>%
    verify(all(str_detect(target, "\\|", negate = TRUE))) %>%
    mutate(otu = glue("{species}|{accession}|{target}"))
  
}

# Remove rogues ----

#' Build a BLAST database.
#'
#' This is a wrapper for makeblastdb.
#'
#' @param in_seqs Character vector of length one; the path to the fasta
#' file containing the sequences to be included in the database.
#' @param db_type Character vector of length one; "nucl" for DNA or "prot"
#' for amino acids (proteins).
#' @param out_name Character vector of length one; name of BLAST database
#' to be created (optional). This will be used to name all database files;
#' if omitted, the name of the `in_seqs` file will be used instead.
#' @param title Character vector of length one; title of BLAST database
#' to be created (optional).
#' @param parse_seqids Logical; should the makeblastdb flag
#' "parse_seqids" be used?
#' @param wd Character vector of length one; working directory. The blast
#' database will be made here.
#' @param echo Logical; should standard error and output be printed?
#' @param ... Additional other arguments. Not used by this function, but
#' meant to be used by \code{\link[drake]{drake_plan}} for tracking
#' during workflows.
#' @return A series of files starting with \code{out_name} and ending in
#' .phr, .pin, .pog, .psd, .psi, psq (for proteins) or .nhr, .nin, .nog,
#' .nsd, .nsi, and .nsq (for DNA) that constitute the BLAST database in
#' the working directory.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK279690/}
#' @examples
#' library(ape)
#' data(woodmouse)
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#' ape::write.FASTA(woodmouse, fs::path(temp_dir, "woodmouse.fasta"))
#' list.files(temp_dir)
#' build_blast_db(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   title = "test db",
#'   out_name = "wood",
#'   parse_seqids = TRUE,
#'   wd = temp_dir)
#' list.files(temp_dir)
#' fs::file_delete(temp_dir)
#' @export
build_blast_db <- function (in_seqs,
                            db_type = "nucl",
                            out_name = NULL,
                            title = NULL,
                            parse_seqids = FALSE,
                            echo = TRUE,
                            wd, ...) {

  # Check input
  assertthat::assert_that(assertthat::is.string(in_seqs))
  assertthat::assert_that(assertthat::is.string(db_type))
  assertthat::assert_that(assertthat::is.string(wd))
  assertthat::assert_that(is.logical(parse_seqids))
  assertthat::assert_that(assertthat::is.string(title) | is.null(title))
  assertthat::assert_that(assertthat::is.string(out_name) | is.null(out_name))

  assertthat::assert_that(db_type %in% c("nucl", "prot"))

  wd <- fs::path_abs(wd)
  assertthat::assert_that(assertthat::is.dir(wd))

  in_seqs <- fs::path_abs(in_seqs)
  assertthat::assert_that(assertthat::is.readable(in_seqs))

  # Prepare arguments
  parse_seqids <- if(isTRUE(parse_seqids)) "-parse_seqids" else NULL
  title <- if(!is.null(title)) c("-title", title) else NULL
  out_name <- if(!is.null(out_name)) c("-out", out_name) else NULL

  arguments <- c("-in", in_seqs,
                 "-dbtype", db_type,
                 parse_seqids,
                 out_name,
                 title)

  # run command
  processx::run("makeblastdb", arguments, wd = wd, echo = echo)

}

#' Run a blastn query.
#'
#' This is a wrapper for blastn.
#'
#' @param query Character vector of length one; the path to the fasta
#' file to use as the query sequence(s).
#' @param database Character vector of length one; the name of the blast
#' database.
#' @param out_file Character vector of length one; the name to use for
#' the results file.
#' @param outfmt Character vector of length one; value to pass to
#' \code{blastn} \code{outfmt} argument. Default = "6".
#' @param other_args Character vector; other arguments to pass on to
#' \code{blastn}.
#' Must be formatted so that each argument name and its value are
#' separate, consecutive elements of the vector, e.g.,
#' \code{c("-evalue", 10, "-num_threads", 1)}.
#' The argument name must be preceded by a hyphen.
#' For a list of options, run \code{blastn -help}.
#' @param wd Character vector of length one; working directory. The blast
#' search will be conducted here.
#' @param echo Logical; should standard error and output be printed?
#' @param ... Additional other arguments. Not used by this function,
#' but meant to be used by \code{\link[drake]{drake_plan}} for tracking
#' during workflows.
#' @return A tab-separated text file with the results of the blastn
#' search, named with the value of \code{out_file}.
#' @author Joel H Nitta, \email{joelnitta@@gmail.com}
#' @references \url{https://www.ncbi.nlm.nih.gov/books/NBK279690/}
#' @examples
#' library(ape)
#'
#' # Make temp dir for storing files
#' temp_dir <- fs::dir_create(fs::path(tempdir(), "baitfindR_example"))
#'
#' # Write out ape::woodmouse dataset as DNA
#' data(woodmouse)
#' ape::write.FASTA(woodmouse, fs::path(temp_dir, "woodmouse.fasta"))
#'
#' # Make blast database
#' build_blast_db(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   db_type = "nucl",
#'   out_name = "wood",
#'   parse_seqids = TRUE,
#'   wd = temp_dir)
#'
#' # Blast the original sequences against the database
#' blast_n(
#'   fs::path(temp_dir, "woodmouse.fasta"),
#'   database = "wood",
#'   out_file = "blastn_results",
#'   wd = temp_dir,
#'   echo = TRUE
#' )
#'
#' # Take a look at the results.
#' readr::read_tsv(
#'   fs::path(temp_dir, "blastn_results"),
#'   col_names = FALSE
#'   )
#'
#' # Cleanup.
#' fs::file_delete(temp_dir)
#' @export
blast_n <- function (query,
                     database,
                     out_file = NULL,
                     outfmt = "6",
                     other_args = NULL,
                     echo = TRUE,
                     wd,
                     ...) {

  # Check input

  assertthat::assert_that(assertthat::is.string(query))
  assertthat::assert_that(assertthat::is.string(database))
  assertthat::assert_that(assertthat::is.string(out_file) | is.null(out_file))
  assertthat::assert_that(assertthat::is.string(outfmt))
  assertthat::assert_that(is.character(other_args) | is.null(other_args))
  assertthat::assert_that(is.logical(echo))
  assertthat::assert_that(assertthat::is.string(wd))
  assertthat::assert_that(
    length(other_args) > 1 | is.null(other_args),
    msg = "other_args not formatted correctly.")

  wd <- fs::path_abs(wd)
  assertthat::assert_that(assertthat::is.dir(wd))

  query <- fs::path_abs(query)
  assertthat::assert_that(assertthat::is.readable(query))

  # modify arguments
  if(!is.null(out_file)) out_file <- c("-out", out_file)

  arguments <- c("-query", query,
                 "-db", database,
                 "-outfmt", outfmt,
                 out_file,
                 other_args)

  # run command
  processx::run("blastn", arguments, wd = wd, echo = echo)

}

#' Convert a list-column of DNA sequences in a tibble to a list of class DNAbin
#'
#' @param seqtbl Tibble containing one DNA sequence per row
#' in a list-column
#' @param name_col Name of column with sequence name
#' @param seq_col Name of column with sequences (list-column)
#' @return List of class DNAbin
#' 
seqtbl_to_dnabin <- function(seqtbl, name_col = "accession", seq_col = "seq") {
  require(ape)
  
  # Check input
  assertthat::assert_that(inherits(seqtbl, "tbl"))
  assertthat::assert_that(assertthat::is.string(name_col))
  assertthat::assert_that(name_col %in% colnames(seqtbl))

  # Extract sequences from metadata and rename
  seqs_dnabin <- do.call(c, seqtbl[[seq_col]])
  names(seqs_dnabin) <- seqtbl[[name_col]]
  
  # Make sure that went OK
  assertthat::assert_that(is.list(seqs_dnabin))
  assertthat::assert_that(inherits(seqs_dnabin, "DNAbin"))
  assertthat::assert_that(all(names(seqs_dnabin) == seqtbl[[name_col]]))

  seqs_dnabin
}

#' Convert DNA sequences from list of class DNAbin to tibble
#'
#' @param dnabin List or matrix of class DNAbin
#' @param name_col Name of column to use for sequence name
#' @param seq_col Name of column to use for sequences (list-column)
#' @return Tibble with one row per sequence
#' 
dnabin_to_seqtbl <- function(dnabin, name_col = "accession", seq_col = "seq") {
  # Convert to a list in case DNA sequences are aligned (in matrix)
  dnabin <- as.list(dnabin)

  # Check input
  assertthat::assert_that(inherits(dnabin, "DNAbin"))
  assertthat::assert_that(assertthat::is.string(name_col))
  assertthat::assert_that(assertthat::is.string(seq_col))

  tibble::tibble(
    "{seq_col}" := split(dnabin, 1:length(dnabin)), 
    "{name_col}" := names(dnabin))
}

#' Make a blast database for ferns
#'
#' @param seqtbl Tibble; fern sequences with column `otu` and `seq`
#' @param blast_db_dir Folder to write BLAST database
#' @param out_name Name of BLAST database
#'
#' @return Paths to components of BLAST database. Externally, database will be created
#' 
make_fern_blast_db <- function(seqtbl, blast_db_dir, out_name) {

	# Extract sequences from metadata and rename
	fern_seqs <- seqtbl_to_dnabin(seqtbl, "otu")
	
	# Remove any gaps
	fern_seqs <- ape::del.gaps(fern_seqs)
	
	# Write out sequences to temporary file
	fern_seqs_path <- tempfile(
		pattern = digest::digest(fern_seqs),
		fileext = ".fasta") %>%
		fs::path_abs()
	if(fs::file_exists(fern_seqs_path)) {fs::file_delete(fern_seqs_path)}
	ape::write.FASTA(fern_seqs, fern_seqs_path)

  # Define output files
	out_files <- fs::path(blast_db_dir, out_name) %>%
		paste0(c(".nhr", ".nin", ".nsq"))

  # Delete any existing output
  for (i in seq_along(out_files)) {
    if(fs::file_exists(out_files[[i]])) {
      fs::file_delete(out_files[[i]])
    }
  }
		
	# Create blast DB (side-effect)
	build_blast_db(                                                     
		fern_seqs_path,                            
		title = out_name,                                                
		out_name = out_name,                                                
		parse_seqids = FALSE,                                              
		wd = blast_db_dir)
	
  # Return path to output file
  out_files

}

#' Run BLAST on DNA sequences in a tibble
#'
#' @param seqtbl Tibble containing DNA sequences as a list-column
#' @param name_col Column with sequence names
#' @param seq_col Column with DNA sequences
#' @param blastdb_files Full paths to BLAST database files
#' @param max_target_seqs Number of maximum hits to return per query
#'
#' @return Tibble
#' 
blast_seqtbl <- function (
  seqtbl, name_col = "otu", seq_col = "seq", 
  blastdb_files, max_target_seqs = 10) {

  # Convert sequences to DNAbin
  query_seqs <- seqtbl_to_dnabin(seqtbl, name_col = name_col, seq_col = seq_col)

  # Remove any gaps
	query_seqs <- ape::del.gaps(query_seqs)
	
	# Write out sequences to temporary file
	query_seqs_path <- tempfile(
		pattern = digest::digest(query_seqs),
		fileext = ".fasta") %>%
		fs::path_abs()
	if(fs::file_exists(query_seqs_path)) {fs::file_delete(query_seqs_path)}
	ape::write.FASTA(query_seqs, query_seqs_path)

  # Define name of output file
  blast_out_path <- tempfile(
		pattern = digest::digest(query_seqs),
		fileext = ".csv") %>%
		fs::path_abs()
	if(fs::file_exists(blast_out_path)) {fs::file_delete(blast_out_path)}

  # Query sequences
  blast_n(
    query = query_seqs_path,
    database = fs::path_file(blastdb_files) %>% fs::path_ext_remove() %>% unique(),
    out_file = blast_out_path, # output as tsv format '6'
    other_args = c("-max_target_seqs", max_target_seqs),
    wd = fs::path_dir(blastdb_files) %>% unique(),
    echo = TRUE
  )
  
  # Read in sequences.
  # BLAST doesn't output column headers, so we need to specify 
  # (make sure they match correctly first!)
  fmt6_cols <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  
  # Read in BLAST output
  blast_results <- readr::read_tsv(
    blast_out_path,
    col_names = fmt6_cols,
    col_types = "ccdddddddddd" # two ID cols are char, rest is numeric
  )
  
  # Cleanup
  if(fs::file_exists(blast_out_path)) fs::file_delete(blast_out_path)
  if(fs::file_exists(query_seqs_path)) fs::file_delete(query_seqs_path)
  
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
#' 
#' @return Tibble; list of rogue sequences showing which families were matched
detect_rogues <- function(metadata_with_seqs, blast_results, ppgi) {
  
  ### Detect rogues ###
  
  ### Check for monotypic families within each gene ###
  # These by default will match to a different (non-self)
  # family, so are not valid for checking rogues.
  exclude_from_rogues <-
    metadata_with_seqs %>% 
    dplyr::left_join(
      select(ppgi, genus, family), 
      by = "genus") %>%
    assert(not_na, genus, family) %>%
    dplyr::add_count(family, target) %>%
    dplyr::filter(n == 1) %>%
    dplyr::mutate(
      otu = glue("{species}|{accession}|{target}") %>% stringr::str_replace_all(" ", "_")
    ) %>%
    select(otu, family, target)
  
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
    verify(all(str_count(qseqid, "\\|") == 2)) %>%
    mutate(
      accession = str_split(qseqid, "\\|") %>% map_chr(2),
      target = str_split(qseqid, "\\|") %>% map_chr(3)) 
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

#' Parse voucher specimen data in data downloaded from GenBank
#' 
#' @param genbank_seqs_tibble Tibble; metadata for GenBank sequences
#' including columns for `gene`, `species`, `specimen_voucher`, 
#' `accession` (GenBank accession), and `length` (gene length in bp)
#'
#' @return Tibble including column `specimen_voucher`
#' 
parse_voucher <- function(genbank_seqs_tibble) {
  genbank_seqs_tibble %>%
    assert(is_uniq, otu) %>%
    # First, filter to only records with specimen voucher (will add others later)
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
      specimen_voucher_raw = map2(meta_data, meta_sel_index, ~magrittr::extract(.x, .y))
    ) %>%
    unnest(specimen_voucher_raw) %>%
    select(subtype, subname, specimen_voucher_raw) %>%
    # Join back onto original data (so, adds `specimen_voucher_raw` column)
    right_join(genbank_seqs_tibble, by = c("subtype", "subname")) %>%
    select(target, species, specimen_voucher_raw, publication, accession, seq_len, otu) %>%
    # In some cases, there are multiple vouchers including `s.n.` and 
    # a numbered voucher. Use only the numbered voucher.
    # (still have some cases with multiple vouchers per accession though
    #  eg., MH101453, KY711736)
    assert(not_na, target, species, accession) %>%
    add_count(target, species, accession) %>%
    filter(!(n > 1 & str_detect(specimen_voucher_raw, "s\\.n\\."))) %>%
    select(-n) %>%
    # Modify voucher:
    # Drop part of voucher name in parentheses
    # this is typically the herbarium; we just want to match on collection,
    # not herbarium where it is stored
    # Also drop "et al"
    mutate(
      specimen_voucher = 
        str_remove_all(specimen_voucher_raw, "\\([^\\(|^\\)]+\\)") %>%
        str_remove_all("et al\\.") %>%
        str_remove_all("et al") %>%
        # Drop characters that are not alphanumeric
        str_remove_all("[^a-zA-Z0-9]+") %>%
        # Convert to lowercase
        str_to_lower() %>%
        str_squish() %>%
        # Replace empty string with NA
        dplyr::na_if("")
    ) %>%
    # Check that all OTUs are accounted for
    # (some are duplicated if there are multiple vouchers)
    verify(all(otu %in% genbank_seqs_tibble$otu)) %>%
    verify(all(genbank_seqs_tibble$otu %in% otu))
}

#' Select a set of GenBank genes by joining on various criteria
#'
#' It would be preferable to join accessions across genes based on specimen
#' voucher, but most accessions lack voucher data.
#' So joining is done based on any of three criteria:
#' - 1. Join across "monophyletic" species (confirmed monophyletic within each
#'      target locus)
#' - 2. Join by voucher
#' - 3. Join by publication, if only one publication for that species
#' - 4. (don't join if none of these conditions are met)
#'
#' After joining accessions, filter to only one set of sequences per species
#' (species), prioritizing in order:
#' - 1. those with rbcL + other genes
#' - 2. those with rbcL
#' - 3. those with the greatest combined length of other genes
#'
#' @param sanger_seqs_with_voucher_data Tibble; metadata for GenBank sequences
#' including columns for `gene`, `species`, `specimen_voucher`,
#' `accession` (GenBank accession), and `length` (gene length in bp)
#' @param mpcheck_monophy Results of species monophyly check
#
#' @return Tibble in wide format joining genes based on species + voucher
#'
select_genbank_genes <- function(sanger_seqs_with_voucher_data,
  mpcheck_monophy) {

  # Check that input names match arguments
  check_args(match.call())

  ### Some pre-processing ###
  # Check overall monophyly by species:
  # must be monophyletic across all genes
  # for each gene with >1 sequence
  mpcheck_monophy_overall <- mpcheck_monophy %>%
    filter(!is.na(is_monophy)) %>%
    select(species, is_monophy) %>%
    group_by(species) %>%
    summarize(is_monophy = all(is_monophy)) %>%
    ungroup()

  # Add column with monophyly data
  genbank_seqs_tibble_with_specimen_dat <- #nolint
    sanger_seqs_with_voucher_data %>%
    left_join(mpcheck_monophy_overall, by = "species") %>%
    # Convert seq len to integer for comparisons
    mutate(seq_len = as.integer(seq_len)) %>%
    assert(not_na, seq_len) %>%
    assert(within_bounds(0, Inf, include.lower = TRUE), seq_len)

  ### Join across accessions ###
  # Don't select any species yet, just do joins.
  # So the species will overlap between join sets.

  # Join monophyletic taxa based on species name
  genbank_seqs_tibble_wide_monophy <- #nolint
    # Filter to only monophyletic taxa
    genbank_seqs_tibble_with_specimen_dat %>%
    filter(is_monophy == TRUE) %>%
    # Make sure we have at least one accession
    verify(
      nrow(.) > 1,
      error_fun = err_msg("No monophyletic species detected")) %>%
    assert(not_na, seq_len, is_monophy, species, target) %>%
    # Select single longest sequence per species per target
    group_by(species, target) %>%
    slice_max(order_by = seq_len, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # Convert to wide format, joining on species
    # - check that joining conditions are unique
    assert_rows(col_concat, is_uniq, target, species) %>%
    # - first split into a list of dataframes by target
    group_by(target) %>%
    group_split() %>%
    # For each dataframe, convert to wide format and
    # rename "seq_len" and "accession" columns by target
    map(
      ~pivot_wider(.,
                   id_cols = c("species",),
                   names_from = target,
                   values_from = c("seq_len", "accession")
      )
    ) %>%
    # Join the target sequences by species
    reduce(full_join, by = "species", na_matches = "never") %>%
    mutate_at(vars(contains("seq_len")), ~replace_na(., 0)) %>%
    # Add total length of all targets (need to sum row-wise)
    # https://stackoverflow.com/questions/31193101/how-to-do-rowwise-summation-over-selected-columns-using-column-index-with-dplyr #nolint
    mutate(total_seq_len = pmap_dbl(select(., contains("seq_len")), sum)) %>%
    mutate(join_by = "monophy")

  # Join based on voucher
  genbank_seqs_tibble_wide_voucher <- #nolint
    genbank_seqs_tibble_with_specimen_dat %>%
    # Exclude specimens lacking a voucher
    filter(!is.na(specimen_voucher)) %>%
    # Make sure we have at least one accession
    verify(
      nrow(.) > 1,
      error_fun = err_msg("No species detected with a voucher")) %>%
    assert(not_na, seq_len, target, specimen_voucher, species) %>%
    # Select single longest sequence per voucher per species per target
    group_by(target, specimen_voucher, species) %>%
    slice_max(order_by = seq_len, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # Exclude if there is only one target gene for the voucher
    # (nothing to join)
    add_count(specimen_voucher, species) %>%
    filter(n > 1) %>%
    select(-n) %>%
    # Convert to wide format
    assert_rows(col_concat, is_uniq, target, specimen_voucher, species) %>%
    group_by(target) %>%
    group_split() %>%
    map(
      ~pivot_wider(.,
                   id_cols = c("species", "specimen_voucher"),
                   names_from = target,
                   values_from = c("seq_len", "accession")
      )
    ) %>%
    # Join the target sequences by species + voucher
    reduce(
      full_join, by = c("species", "specimen_voucher"),
      na_matches = "never") %>%
    mutate_at(vars(contains("seq_len")), ~replace_na(., 0)) %>%
    # Add total length of all targets
    mutate(total_seq_len = pmap_dbl(select(., contains("seq_len")), sum)) %>%
    mutate(join_by = "voucher")

  # Join based on publication conditions: must have only one publication per
  # species across target loci
  genbank_seqs_tibble_wide_publication <- #nolint
    genbank_seqs_tibble_with_specimen_dat %>%
    # Filter to only those with publication data
    filter(!is.na(publication)) %>%
    # Make sure we have at least one accession
    verify(
      nrow(.) > 1,
      error_fun = err_msg("No species detected with a publication")) %>%
    # Convert publication to lower case to account for
    # mistakes in capitalization
    mutate(publication = str_to_lower(publication)) %>%
    # Filter to those with one publication per species
    assert(not_na, species, target, publication) %>%
    add_count(species, target) %>%
    filter(n == 1) %>%
    select(-n) %>%
    # Of these, filter to only taxa with multiple target loci
    add_count(species) %>%
    filter(n > 1) %>%
    select(-n) %>%
    # Convert to wide format
    assert_rows(col_concat, is_uniq, target, species, publication) %>%
    group_by(target) %>%
    group_split() %>%
    map(
      ~pivot_wider(.,
                   id_cols = c("species", "publication"),
                   names_from = target,
                   values_from = c("seq_len", "accession")
      )
    ) %>%
    # Join the target sequences by species + publication
    reduce(full_join,
      by = c("species", "publication"),
      na_matches = "never") %>%
    mutate_at(vars(contains("seq_len")), ~replace_na(., 0)) %>%
    # Add total length of all targets
    mutate(total_seq_len = pmap_dbl(select(., contains("seq_len")), sum)) %>%
    mutate(join_by = "publication")

  # Make vector of non-rbcL accession columns for selection
  # (can't use negative look-ahead in select() and friends)
  # https://github.com/r-lib/tidyselect/issues/58
  other_acc_cols <- colnames(genbank_seqs_tibble_wide_monophy) %>%
    magrittr::extract(
      str_detect(., "accession_(?!.*rbcL)")
    )

  # Combine wide joined tibbles
  joined_all <- bind_rows(
    genbank_seqs_tibble_wide_monophy,
    genbank_seqs_tibble_wide_voucher,
    genbank_seqs_tibble_wide_publication
  ) %>%
  rowwise() %>%
  mutate(
    # Add column indicating presence of any gene other than rbcL
    accession_other = any(!is.na(c_across(any_of(other_acc_cols))))
    ) %>%
  ungroup() %>%
  mutate(
    # Add column indicating presence of rbcL and another gene
    rbcL_and_other = case_when(
      !is.na(accession_rbcL) & accession_other == TRUE ~ TRUE,
      TRUE ~ FALSE
    )
  )

  ### Final sequence selection ###

  # 1. Filter to best joined sequences per species with rbcL and another gene
  rbcl_joined <- joined_all %>%
    filter(rbcL_and_other == TRUE) %>%
    # Slice to longest total seq per species
    group_by(species) %>%
    slice_max(n = 1, order_by = total_seq_len, with_ties = FALSE) %>%
    ungroup()

  # 2. Of the remainder, take the best unjoined rbcL sequence per species
  rbcl_unjoined <-
    # Start with all seqs in long format
    genbank_seqs_tibble_with_specimen_dat %>%
    # Filter to rbcL
    filter(target == "rbcL", seq_len > 0, !is.na(seq_len)) %>%
    # Exclude species with joined rbcL
    anti_join(rbcl_joined, by = "species") %>%
    # Slice to longest rbcL per species
    group_by(species) %>%
    slice_max(n = 1, order_by = seq_len, with_ties = FALSE) %>%
    ungroup() %>%
    # Make sure species are unique
    assert(is_uniq, species) %>%
    # Format like other wide data
    rename(accession_rbcL = accession, seq_len_rbcL = seq_len) %>%
    mutate(total_seq_len = seq_len_rbcL) %>%
    mutate(join_by = "unjoined") %>%
    select(-target, -specimen_voucher_raw, -otu, -is_monophy)

  # Prepare to select remainder
  # - Select best species for "other" joined genes
  other_joined <-
    joined_all %>%
    # Exclude species that have rbcL
    anti_join(rbcl_unjoined, by = "species") %>%
    anti_join(rbcl_joined, by = "species") %>%
    # Slice to longest total seq per species
    group_by(species) %>%
    slice_max(n = 1, order_by = total_seq_len, with_ties = FALSE) %>%
    ungroup()

  # - Select best species for "other" unjoined genes
  # (will overlap in species with other_joined)
  other_unjoined <-
    genbank_seqs_tibble_with_specimen_dat %>%
    # Exclude species that have rbcL
    anti_join(rbcl_unjoined, by = "species") %>%
    anti_join(rbcl_joined, by = "species") %>%
    # Filter to single longest combine sequence length
    group_by(species) %>%
    slice_max(n = 1, order_by = seq_len, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(join_by = "unjoined") %>%
    select(-specimen_voucher_raw, -otu, -is_monophy) %>%
    # Convert to wide format, joining on species
    # - check that joining conditions are unique
    assert_rows(col_concat, is_uniq, target, species) %>%
    # - first split into a list of dataframes by target
    group_by(target) %>%
    group_split() %>%
    # For each dataframe, convert to wide format and
    # rename "seq_len" and "accession" columns by target
    map(
      ~pivot_wider(.,
                   id_cols = c("species",),
                   names_from = target,
                   values_from = c("seq_len", "accession")
      )
    ) %>%
    # Join the target sequences by species
    reduce(full_join, by = "species", na_matches = "never") %>%
    mutate_at(vars(contains("seq_len")), ~replace_na(., 0)) %>%
    mutate(join_by = "unjoined") %>%
    # Add total length of all targets (need to sum row-wise)
    mutate(total_seq_len = pmap_dbl(select(., contains("seq_len")), sum)) %>%
    assert(is_uniq, species)

  # 3. Make final selection by total sequence length from unjoined and joined
  # "other" genes
  # (since the joined sequence may be shorter than an unjoined sequence for a
  # given species: e.g., sp1 unjoined atpA 1200 bp vs. sp1 joined trnlF + rps4 500 bp) #nolint
  other_genes <-
  bind_rows(
    other_unjoined,
    other_joined
  ) %>%
    assert(not_na, total_seq_len) %>%
    assert(within_bounds(0, Inf, include.lower = TRUE), total_seq_len) %>%
    # Filter to single longest combine sequence length
    group_by(species) %>%
    slice_max(n = 1, order_by = total_seq_len, with_ties = FALSE) %>%
    ungroup()

  # Combine final selected sequences
  rbcl_joined %>%
    bind_rows(rbcl_unjoined) %>%
    bind_rows(other_genes) %>%
    select(-rbcL_and_other, -accession_other) %>%
    # Make sure species are unique
    assert(not_na, species, join_by) %>%
    assert(is_uniq, species) %>%
    # Make sure all original input taxa are present
    verify(all(species %in% sanger_seqs_with_voucher_data$species)) %>%
    verify(all(sanger_seqs_with_voucher_data$species %in% species)) %>%
    arrange(species)
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

#' Inspect rogue sequences
#' 
#' Sanger rogue sequences were assessed by checking the top BLAST hit from
#' an all-by-all BLAST. Those whose families don't match were flagged as
#' rogues. However, some of this may be due to incorrect taxonomic treatment.
#' This checks the family of the original GenBank accession to see if the
#' taxonomic treatment may need to change.
#'
#' @param sanger_seqs_rogues Tibble; Potential rogue sequences identified
#' by query family not matching hit family.
#' @param raw_meta_all Tibble; metadata of Sanger sequences including NCBI
#' taxonomic ID (taxid) and GenBank accession number.
#' @param ncbi_names_query Tibble; Names downloaded from NCBI taxonomy database
#' @param ppgi_taxonomy Tibble; PPGI taxonomy.
#'
#' @return Tibble with potential "real rogue" sequences flagged as 1
#'
inspect_rogues <- function(
	sanger_seqs_rogues,
	raw_meta_all,
	ncbi_names_query,
	ppgi_taxonomy) {

  # Check that input names match arguments
  check_args(match.call())
	
	# Filter original GenBank names to accepted name
	gb_names <-
		ncbi_names_query %>%
		filter(accepted == TRUE) %>%
		# Combine `scientific_name`, `species` filling in NAs as needed
		mutate(gb_name = dplyr::coalesce(scientific_name, species)) %>%
		select(taxid, gb_name) %>%
		unique()
	
	sanger_seqs_rogues %>%
		# Split query text into component parts
		separate(qseqid, c("species", "accession", "gene"), sep = "\\|") %>%
		verify(all(gene == target)) %>%
		select(-target) %>%
		# Add taxid
		left_join(
			unique(select(raw_meta_all, accession, taxid)), by = "accession"
		) %>%
		# Join original GenBank names by taxid
		assert(not_na, taxid) %>%
		left_join(gb_names, by = "taxid") %>%
		assert(not_na, gb_name) %>%
		# Join original GenBank family
		mutate(gb_genus = str_split(gb_name, " ") %>%
					 	map_chr(1)) %>%
		left_join(
			select(ppgi_taxonomy, gb_genus = genus, gb_family = family), 
        by = "gb_genus"
		) %>%
		select(-gb_genus) %>%
		# Consider a potential "real rogue" if
		# the query family and genbank family agree
		mutate(
			real_rogue = case_when(
				q_family == gb_family ~ 1,
				TRUE ~ 0
			)
		)
	
}

# Check for species monophyly in Sanger loci ----

#' Check monophyly of a species
#' 
#' Monophyly check done with ape::is.monophyletic(), which is rather
#' slow. Speed things up by running in parallel.
#'
#' @param mpcheck_sliced Tibble (seqtbl) with accession and species names
#' @param mpcheck_tree Phylogenetic tree for a single target locus (or
#' dataframe containing this in column 'tree')
#' @param workers Number of workers to run in parallel
#'
#' @return Tibble, with logical column `is_monophy` indicating monophyly
#' of species with >1 accession. If species has only 1 accession, `is_monophy` 
#' is `NA`.
check_monophy <- function(mpcheck_sliced, mpcheck_tree, workers) {
	
	# Change back to sequential when done (including on failure)
	on.exit(future::plan(future::sequential), add = TRUE)
	
	# Set backend for parallelization
	future::plan(future::multisession, workers = workers)

  # If mpcheck_tree input is dataframe, extract tree
  if(inherits(mpcheck_tree, "data.frame")) {
    assertthat::assert_that(nrow(mpcheck_tree) == 1)
    assertthat::assert_that("tree" %in% colnames(mpcheck_tree))
    mpcheck_tree <- mpcheck_tree$tree[[1]]
  }
	
  # Make dataframe of taxa to check for monophyly (>1 accession)
	species_to_check <-
		mpcheck_sliced %>% 
    verify(all(accession %in% mpcheck_tree$tip.label)) %>%
		add_count(species) %>%
		filter(n > 1) %>%
		select(species, accession) %>%
		unique()
	
  # and dataframe of accessions to skip (1 accession each)
	species_to_skip <- mpcheck_sliced %>% 
		add_count(species, name = "n_accs") %>%
		filter(n_accs == 1) %>%
		select(species, n_accs)
	
  # check monophyly
	res <- species_to_check %>%
		group_by(species) %>%
		add_count(name = "n_accs") %>%
		nest(data = accession) %>%
		ungroup() %>%
		mutate(
      # ape::is.monophyletic is slow, so run in parallel
			is_monophy = furrr::future_map_lgl(
				data, 
				~ape::is.monophyletic(
					phy = mpcheck_tree, 
					tips = .x$accession, reroot = TRUE, plot = FALSE),
        .options = furrr::furrr_options(seed = TRUE)
			)
		) %>%
		select(species, n_accs, is_monophy) %>%
		bind_rows(species_to_skip) %>%
    mutate(target = unique(mpcheck_sliced$target))
	
	# Close parallel workers
	future::plan(future::sequential)
	
	res
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

#' Format a query string to download fern plastome sequences from GenBank
#'
#' There is no strict definition of plastome, so consider anything >7000 bp
#' to be "plastome"
#'
#' @param start_date String; start date to include sequences
#' @param end_date String; end date to include sequences
#'
#' @return String
format_fern_plastome_query <- function(start_date = "1980/01/01", end_date, strict = FALSE) {
  if(isTRUE(strict)) {
    # require "partial or complete" "genome" if strict
    # exclude one accession that causes problems for extracting genes AP004638
    query <- glue('Polypodiopsida[ORGN] AND (plastid OR chloroplast) AND 7000:500000[SLEN] AND ("{start_date}"[PDAT]:"{end_date}"[PDAT]) AND (partial OR complete) AND genome NOT AP004638') # no lint
  } else {
    query <- glue('Polypodiopsida[ORGN] AND (plastid OR chloroplast) AND 7000:500000[SLEN] AND ("{start_date}"[PDAT]:"{end_date}"[PDAT])') # no lint
  }
  query
}

#' Download plastome metadata for ferns
#' 
#' Using `strict = TRUE` will (most likely) restrict search to properly
#' annotated plastomes, but fewer sequences. `strict = FALSE` will result in
#' more hits, but some may not be annotated (at all)
#' 
#' @param start_date Earliest date to download
#' @param end_date Most recent date to download
#' @param strict Logical; use strict version of search string or not?
#' @param outgroups Dataframe of genbank accession numbers
#' and species names to use for outgroups (non-fern taxa that wouldn't
#' be captured by the query)
#' @param accs_exclude Character vector of GenBank accessions to exclude
#'
#' @return Tibble of GenBank metadata combining the results of
#' the query and the outgroups.
#' 
download_plastome_metadata <- function (start_date = "1980/01/01", end_date, 
  strict = FALSE, outgroups, accs_exclude = NULL) {
  
  assertthat::assert_that(assertthat::is.string(start_date))
  assertthat::assert_that(assertthat::is.string(end_date))
  
  # Format GenBank query: all ferns plastomes within specified dates
  # There is no formal category for "whole plastome" in genbank, so use size
  # cutoff: >7000 and < 500000 bp
  ingroup_query = format_fern_plastome_query(start_date = start_date,
    end_date = end_date, strict = strict)

  # Fetch standard metadata
  ingroup_metadata_raw <- fetch_metadata(
    query = ingroup_query,
    col_select = c("gi", "caption", "taxid", "title", "subtype", "subname", "slen"))
  
  ingroup_metadata <-
    ingroup_metadata_raw %>%
    rename(accession = caption) %>%
    # GenBank accession should be non-missing, unique
    assert(not_na, accession) %>%
    assert(is_uniq, accession) %>%
    # Remove GenBank duplicate sequences starting with "NC_"
    # - consider accessions with same species and seq length to be same
    # if they only differ in accession.
    mutate(maybe_dup = case_when(
      str_detect(accession, "NC_") ~ 1,
      TRUE ~ 0)) %>%
    group_by(taxid, slen) %>%
    arrange(maybe_dup) %>%
    slice(1) %>%
    ungroup %>%
    # Sort by accession
    arrange(accession) %>%
    select(-maybe_dup)
  
  # Download outgroup metadata: use one representative of each major 
  # seed plant, lycophyte, and bryo group
  og_query <- outgroups %>% pull(accession) %>% paste(collapse = "[accession] OR ") %>%
    paste("[accession]", collapse = "", sep = "")
  
  outgroup_metadata_raw <- fetch_metadata(
    query = og_query,
    # don't fetch `slen` (length of accession; will calculate length of actual sequence later)
    col_select = c("gi", "caption", "taxid", "title", "subtype", "subname", "slen"))
  
  # Combine ingroup and outgroup data
  outgroup_metadata <-
    outgroup_metadata_raw %>%
    rename(accession = caption) %>%
    verify(all(accession %in% outgroups$accession)) %>%
    verify(all(outgroups$accession %in% accession))
  
  combined_metadata <- bind_rows(ingroup_metadata, outgroup_metadata)

  # Exclude accessions in accs_exclude
  if (!is.null(accs_exclude)) {
    combined_metadata <-
      combined_metadata %>%
      filter(!accession %in% accs_exclude)
  }

  combined_metadata

}

#' Fetch a set of genes from a plastome
#'
#' @param genes Vector of gene names
#' @param accession Genbank accession number (of a plastome)
#' @param max_length Maximum length to accept for genes. Used to filter
#' out any abnormally (probably erroneously) long genes.
#'
#' @return Dataframe with columns for `accession`, `gene`, and `seq`
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
  temp_file <- tempfile(pattern = "gb_records_", fileext = ".txt")
  if(fs::file_exists(temp_file)) fs::file_delete(temp_file)
  
  reutils::efetch(uid, "nucleotide", rettype = "gb", retmode = "text", outfile = temp_file)
  
  # Parse flatfile
  gb_entry <- readr::read_file(temp_file)
  
  # get the results (includes NULL if errored)
  extracted_genes <- map(genes, ~extract_sequence(gb_entry, .)) %>% 
    transpose() %>% 
    magrittr::extract2("seq")
  
  # subset gene names to those to successfully extracted
  genes_successful <- genes[!map_lgl(extracted_genes, is.null)]
  
  # Subset results to successful genes, filter by length, and set names
  extracted_genes_filtered <-
    extracted_genes %>%
    # Drop errors
    compact() %>%
    flatten()
  
  # Make sure length of filtered names matches
  assertthat::assert_that(
    length(extracted_genes_filtered) == length(genes_successful)
  )
  
  extracted_genes_filtered <-
    extracted_genes_filtered %>%
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
  
  # Return as a dataframe, so the results can be combined later
  tibble(
    accession = accession,
    gene = names(extracted_genes_filtered),
    seq = extracted_genes_filtered
    )
}

#' Rename plastome sequences
#'
#' Renames plastome sequences downloaded from GenBank by their resolved names
#' from pteridocat, filters to one best (longest) sequence per genus
#'
#' @param plastome_genes_raw Tibble (seqtbl); gene sequences downloaded from 
#' GenBank
#' @param plastome_metadata_renamed Tibble; metadata of sequences
#' including resolved names.
#'
#' @return Tibble (seqtbl); gene sequences downloaded from 
#' GenBank with column "accession" including species_accession
#'
rename_and_filter_raw_plastome_seqs <- function(plastome_genes_raw, plastome_metadata_renamed) {
  plastome_genes_raw %>%
    # Add resolved names  
    left_join(
      select(plastome_metadata_renamed, species, accession),
      by = "accession"
    ) %>%
    verify(nrow(.) == nrow(plastome_genes_raw)) %>%
    assert(not_na, everything()) %>%
    # Add sequence length
    mutate(
      seqln = map_dbl(seq, ~length(.[[1]])),
      genus = str_split(species, "_") %>% map_chr(1)) %>%
    verify(all(str_detect(genus, " ", negate = TRUE))) %>%
    assert(not_na, everything()) %>%
    group_by(gene, genus) %>%
    # Filter to longest sequence per genus
    slice_max(order_by = "seqln", n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      accession = glue::glue("{species}_{accession}") %>% as.character(),
      seq,
      target = gene
    )
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
#' @param genes_keep Names of genes (loci) to retain regardless of number
#' of accessions recovered
#' @param accs_keep Names of accessions to retain regardless of number of
#' genes recovered
#'
#' @return Dataframe
#'
filter_majority_missing <- function (
  gene_lengths_best,
  genes_keep = c("rps4-trnS", "trnL-trnF"),
  accs_keep) {
  
  # Make sure genes_keep and accs_keep are valid
  assertthat::assert_that(
    all(genes_keep %in% gene_lengths_best$gene),
    msg = "genes_keep not in gene_lengths_best$gene"
  )
  assertthat::assert_that(
    all(accs_keep %in% gene_lengths_best$accession),
    msg = "accs_keep not in gene_lengths_best$accession"
  )
  
  # Count number of accessions recovered per gene and flag accessions to remove
  missing_acc_summ <-
    gene_lengths_best %>%
    # Assume gene is missing if < 10 bp
    mutate(missing = case_when(
      slen < 10 ~ TRUE,
      TRUE ~ FALSE
    )) %>%
    group_by(gene) %>%
    summarize(
      n_missing_accs = sum(missing),
      .groups = "drop"
    ) %>%
    # Set gene to keep if present in >50% of accessions
    # (or in keep list)
    mutate(
      keep_gene = case_when(
        n_missing_accs < 0.5 * n_distinct(gene_lengths_best$accession) ~ TRUE,
        gene %in% genes_keep ~ TRUE,
        TRUE ~ FALSE
      )
    )
  
  # Count number of genes recovered per accession and flag genes to remove
  missing_gene_summ <-
    gene_lengths_best %>%
    mutate(missing = case_when(
      slen < 10 ~ TRUE,
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
        accession %in% accs_keep ~ TRUE,
        TRUE ~ FALSE
      )
    )
  
  # Filter genes
  gene_lengths_best %>%
    left_join(missing_acc_summ, by = "gene") %>%
    left_join(missing_gene_summ, by = "accession") %>%
    filter(keep_gene == TRUE, keep_acc == TRUE) %>%
    select(-keep_gene, -keep_acc) %>%
    assert_rows(col_concat, is_uniq, accession, gene, species)
}

#' Select a list of plastid sequences to use from a list of plastid genes and
#' plastome metadata
#'
#' Filters list to one best accession per species, only keeping genes that are missing
#' < 50% of accessions and accessions missing < 50% of genes'
#'
#' @param plastome_genes_raw Dataframe of plastid genes. Each row is
#' a sequence for a gene for a plastome accession.
#' @param plastome_metadata_renamed Associated plastome metadata (species names)
#' @param fern_plastome_loci_extract_res Output of extract_from_ref_blast(); spacer sequences
#' in plastomes
#'
select_plastome_seqs <- function (plastome_genes_raw, plastome_metadata_renamed, fern_plastome_loci_extract_res) {
  
  # Check that input names match arguments
  check_args(match.call())

  # Extract list of outgroup accessions to keep in gene selection
  outgroup_accs <-
    plastome_metadata_renamed %>%
    filter(outgroup == TRUE) %>%
    assert(not_na, accession) %>%
    pull(accession)

  # Extract sequences from superCRUNCH results
  plastome_genes_raw <-
  fern_plastome_loci_extract_res %>%
    clean_extract_res("dc-megablast") %>%
    left_join(
      select(plastome_metadata_renamed, accession, species),
      by = "accession"
    ) %>%
    # Check no rows are duplicated, all accs have species
    assert(not_na, everything()) %>%
    assert_rows(col_concat, is_uniq, accession, target, species) %>%
    rename(gene = target)
    
  # Make tibble of gene lengths by accession, including species and voucher
  gene_lengths <- 
    plastome_genes_raw %>%
    # Make sure each accession has at least one gene
    add_count(accession, gene) %>%
    verify(all(n > 0)) %>%
    select(-n) %>%
    # Add column for sequence length
    mutate(slen = map_dbl(seq, ~length(.[[1]]))) %>%
    select(-seq) %>%
    assert(not_na, everything())
  
  # Missing genes (length 0) are not in the original sequences list,
  # so add these by crossing all combinations of accession and gene
  gene_lengths <- 
    purrr::cross_df(list(
      gene = gene_lengths$gene %>% unique, 
      accession = gene_lengths$accession %>% unique)) %>%
    left_join(select(gene_lengths, gene, accession, slen), by = c("gene", "accession")) %>%
    left_join(select(gene_lengths, accession, species) %>% unique, by = "accession") %>%
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
    # first get total rel length for each accession, keeping species column
    group_by(accession, species) %>%
    summarize(
      total_rel_len = sum(rel_len),
      .groups = "drop"
    ) %>%
    # then sort by species and keep the one with the greatest length
    group_by(species) %>%
    slice_max(n = 1, order_by = total_rel_len, with_ties = FALSE) %>%
    ungroup
  
  # Assemble final selection: filtered plastome loci, one accession per species
    gene_lengths %>%
      # Make a table of (relative) gene lengths for
      # the best accession per species
      inner_join(select(best_accessions_by_species, accession), by = "accession") %>%
      left_join(max_lengths, by = "gene") %>%
      mutate(rel_len = slen / max_length) %>%
      assert(not_na, everything()) %>%
      assert_rows(col_concat, is_uniq, accession, gene, species) %>%
      # Filter out accessions missing > 50% of genes
      # and genes absent from > 50% of sequences
      filter_majority_missing(accs_keep = outgroup_accs) %>%
      filter(slen > 1) %>%
      # Add sequence
      select(species, accession, gene) %>%
      left_join(
        select(plastome_genes_raw, accession, gene, seq),
        by = c("accession", "gene")) %>%
      assert(not_na, everything()) %>%
      rename(target = gene)

}

#' Trim aligned plastome genes
#'
#'
#' @param plastome_genes_aligned Tibble (seqtbl); aligned plastome genes.
#'
#' @return Trimmed, aligned plastome genes
#' 
trim_plastome_genes <- function(plastome_genes_aligned) {
  plastome_genes_aligned %>% 
     select(seq, target, accession) %>%
     group_by(target) %>%
     nest(data = c(seq, accession)) %>%
     mutate(
       align_trimmed = map(
         data, trimal, 
         other_args = c("-gt", "0.05"), 
         return_seqtbl = FALSE, 
         name_col_in = "accession"
       )) %>%
     ungroup() %>%
     select(-data)
}

# Combine and align Sanger and plastome sequences ----

#' Assign taxonomic clusters to sequences
#'
#' Only meant for spacer regions. Will drop sequences
#' if fewer than four sequences in a cluster
#'
#' @param plastid_genes_unaligned_combined Combined Sanger and plastome sequences
#' @param ppgi_taxonomy PPGI taxonomy
#' @param plastome_metadata_renamed Plastome metadata with resolved taxonomic names
#' @param target_select Name of selected target region
#'
#' @return seqtbl
#'
assign_tax_clusters <- function(
  sanger_accessions_selection,
  sanger_seqs_combined_filtered,
  plastome_seqs_combined_filtered, 
  ppgi_taxonomy, plastome_metadata_renamed,
  target_select) {

  # Combine sanger and plastome sequences
  combine_sanger_plastome(
      sanger_accessions_selection,
      sanger_seqs_combined_filtered,
      plastome_seqs_combined_filtered) %>%
    # Extract sequences for selected target
    filter(target == target_select) %>%
    # Exclude outgroups
    anti_join(
      filter(plastome_metadata_renamed, outgroup == TRUE), by = "species"
    ) %>%
    # Add taxonomy
    mutate(genus = str_split(species, "_") %>% map_chr(1)) %>%
    assert(not_na, genus) %>%
    left_join(
      select(ppgi_taxonomy, genus, family, suborder, order), by = "genus"
    ) %>%
    assert(not_na, family, order) %>%
    # Count number of taxa per family
    add_count(family) %>%
    # Assign cluster as family if >=2 taxa, none if 1
    mutate(
      cluster = case_when(
        n >= 2 ~ family,
        n < 2 ~ "none"
      )
    ) %>%
    select(-c(n, genus, family, suborder, order)) %>%
    assert(not_na, cluster, seq) %>%
    assert(not_null, seq)
}

#' Combine Sanger and plastome sequences, rename sequences by species
#'
#' @param sanger_accessions_selection Dataframe; selection of Sanger accessions in 
#' wide format. Output of select_genbank_genes()
#' @param sanger_seqs_combined_filtered Dataframe; filtered Sanger sequences
#' Output of combine_and_filter_sanger()
#' @param plastome_seqs_combined_filtered Dataframe; filtered plastome sequences.
#' Output of select_plastome_seqs()
#'
#' @return Dataframe; combined Sanger and plastome sequences in long format
#' 
combine_sanger_plastome <- function(
  sanger_accessions_selection,
  sanger_seqs_combined_filtered,
  plastome_seqs_combined_filtered
) {
  
  sanger_accessions_selection %>%
    # Remove any accessions in plastome data
    anti_join(plastome_seqs_combined_filtered, by = "species") %>%
    # Convert to long form
    select(species, matches("accession")) %>%
    pivot_longer(names_to = "target", values_to = "accession", -species) %>%
    mutate(target = str_remove_all(target, "accession_")) %>%
    filter(!is.na(accession)) %>%
    # Add DNA sequences
    left_join(select(sanger_seqs_combined_filtered, accession, seq, target), by = c("accession", "target")) %>%
    # Add plastome sequences
    bind_rows(plastome_seqs_combined_filtered) %>%
    assert(not_na, everything()) %>%
    assert_rows(col_concat, is_uniq, species, target, accession)
    
}

#' Align sequences in a tibble
#'
#' @param seqtbl Tibble containing one DNA sequence per row
#' in a list-column
#' @param name_col Name of column with sequence name
#' @param seq_col Name of column with sequences (list-column)
#'
#' @return Dataframe with aligned sequences. New column logical column "reversed" will
#' be appended indicating if sequence was reversed when aligning or not.
#' 
align_seqs_tbl <- function(seqs_tbl, name_col = "accession", seq_col = "seq") {
  # Extract sequences, convert to ape DNAbin list
  seqs <-
    seqs_tbl %>%
    # Name column must be unique values
    assert(is_uniq, all_of(name_col)) %>%
    seqtbl_to_dnabin(name_col = name_col, seq_col = seq_col)

  # Set aside metadata
  seqs_data <-
    select(seqs_tbl, -all_of(seq_col)) %>%
    # Metadata can't already contain column name "reversed"
    verify(!"reversed" %in% colnames(.))

  # Align sequences
  alignment <- ips::mafft(
    x = seqs,
    options = "--adjustdirection",
    exec = "/usr/bin/mafft")

  # Join metadata back to aligned sequences
  # - make tibble of which seqs were reversed
  reversed_seqs_tbl <-
  alignment %>%
    dnabin_to_seqtbl(name_col = name_col, seq_col = seq_col) %>%
    select(.data[[name_col]]) %>%
    mutate(
      reversed = str_detect(.data[[name_col]], "_R_"),
      !!name_col := str_remove_all(.data[[name_col]], "_R_")
    )

  # - fix names if mafft changed them
  rownames(alignment) <- str_remove_all(rownames(alignment), "_R_")

  # Join metadata back to aligned sequences
  alignment %>%
    dnabin_to_seqtbl(name_col = name_col, seq_col = seq_col) %>%
    left_join(reversed_seqs_tbl, by = name_col) %>%
    left_join(seqs_data, by = name_col)

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
#' @param seqs DNA sequence alignment; matrix of class 'DNAbin' or seqtbl
#' @param echo Optional; should the output of trimal be printed to the
#' screen?
#' @param other_args; Character vector of additional arguments to pass
#' to `trimal`. Must have one element per word that would normally be
# separated by a space on the command line. E.g., `c("-gt", "0.05")`, etc.
#' @param echo Logical; should messages from trimal be printed to screen?
#' @param name_col_in Name of DNA sequence name column if seqs is seqtbl
#' @param seq_col_in Name of DNA sequence column if seqs is seqtbl
#' @param name_col_out Name of DNA sequence name column if output type is seqtbl
#' @param seq_col_out Name of DNA sequence column if output type is seqtbl
#' @param return_seqtbl Logical; should the results be returned as a DNA-sequence tibble (seqtbl)?
#' If FALSE, will return as list of class 'DNAbin'. Default to TRUE if `seqs` is seqtbl.
#'
#' @return DNA sequence alignment with gappy sites removed by trimal
#' 
trimal <- function (
    seqs, other_args = NULL, echo = FALSE, 
    name_col_in = "accession", seq_col_in = "seq", 
    name_col_out = "accession", seq_col_out = "seq", 
    return_seqtbl = inherits(seqs, "tbl")) {

  # Evaluate return type
  return_seqtbl <- return_seqtbl

  # Convert to DNAbin if in tibble
  if(inherits(seqs, "tbl")) {
    seqs <- seqtbl_to_dnabin(seqs, name_col = name_col_in, seq_col = seq_col_in)
  }
  
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
  args <- c("-in", in_file_name, "-out", out_file_name, "-fasta", other_args)
  
  # Run trimal
  results <- processx::run("trimal", args, wd = temp_wd, 
                           echo = echo)
  
  # Read in results
  results <- ape::read.FASTA(fs::path(temp_wd, out_file_name)) %>% as.matrix
  
  # Clean up
  fs::file_delete(fs::path(temp_wd, out_file_name))
  fs::file_delete(fs::path(temp_wd, in_file_name))
  
  # All done
  if(isTRUE(return_seqtbl)) {
    results <- dnabin_to_seqtbl(results, name_col = name_col_out, seq_col = seq_col_out)
  }

  results
  
}

#' Trim spacer regions
#'
#' @param plastid_spacers_aligned_clusters seqtbl of aligned platid spacers by cluster
#'
#' @return Tibble with list-column of trimmed spacers; each row is
#' one cluster per target region
#' 
trim_spacers_by_cluster <- function(plastid_spacers_aligned_clusters) {

  # Check that input names match arguments
  check_args(match.call())

  plastid_spacers_aligned_clusters %>% 
    select(seq, species, target, cluster) %>%
    group_by(target, cluster) %>%
    nest(data = c(seq, species)) %>%
    # trim spacers lightly to keep most gaps
    mutate(
      align_trimmed = map(
        data, trimal, 
        other_args = c("-gt", "0.01"), 
        return_seqtbl = FALSE,
        name_col_in = "species"
      )) %>%
    ungroup() %>%
    select(-data)
}

#' Trim genes (or spacers)
#'
#' @param plastid_aligned seqtbl of aligned platid genes or spacers
#' @param name_col_in Name of column in plastid_aligned to use as sequence
#'   names when trimming with trimal.
#'
#' @return Tibble with list-column of trimmed genes (or spacers); each row is
#' a target gene (or spacer)
#' 
trim_genes <- function(plastid_aligned, name_col_in = "species") {

  plastid_aligned %>% 
    select(seq, species, target, accession) %>%
    group_by(target) %>%
    nest(data = c(seq, species, accession)) %>%
    # trim genes (0.05) slightly more aggressively than spacers (0.01)
    mutate(
      trim_strength = if_else(
        # spacers have hyphen, genes don't
        str_detect(target, "-"), "0.01", "0.05"),
      align_trimmed = map(
        data, trimal, 
        other_args = c("-gt", trim_strength), 
        return_seqtbl = FALSE, 
        name_col_in = name_col_in
      )) %>%
    ungroup() %>%
    select(-data, -trim_strength)
}

#' Align one representative sequence per spacer region
#'
#' @param plastid_spacers_aligned_trimmed Trimmed aligned spacer sequences in a tibble,
#' with alignments in list-column "align_trimmed"
#' @param plastid_spacers_unaligned seqtbl with taxonomic clusters assigned to spacer regions
#' @param target_select Name of target (e.g., "trnL-trnF") to align.
#' @param exclude_terms Grep string. Any clusters with names matching this will be excluded from results.
#'
#' @return Tibble with reverse-complemented sequences. Column "align-trimmed" includes
#' alignments of clusters; column "seq" includes sequences of individual singletons.
#'
align_rep_spacers <- function(plastid_spacers_aligned_trimmed, plastid_spacers_unaligned, target_select, exclude_terms = NULL) {

# Start with aligned spacer seqs, one alignment per row
  plastid_spacers_aligned_trimmed %>%
  # Unnest alignments to sequences (so each row is now a sequence)
  assert(not_na, align_trimmed) %>%
  mutate(align_trimmed = map(align_trimmed, as.list)) %>%
  mutate(align_trimmed = map(align_trimmed, ~split(., names(.)))) %>%
  unnest(align_trimmed) %>%
  rename(seq = align_trimmed) %>%
  mutate(species = map_chr(seq, names)) %>%
  # Add singletons
  bind_rows(filter(plastid_spacers_unaligned, cluster == "none")) %>%
  # Make sure we have filtered to a single target
  filter(target == target_select) %>%
  verify(n_distinct(target) == 1) %>%
  assert(not_na, target) %>%
  # If cluster is "none", group at species level
  assert(not_na, cluster, target, seq, species) %>%
  mutate(cluster_grp = case_when(
    cluster == "none" ~ species,
    TRUE ~ cluster
  )) %>%
  # Remove gaps
  mutate(seq = map(seq, ape::del.gaps)) %>%
  # Filter to longest sequence per cluster
  assert(not_null, seq) %>%
  mutate(seq_len = map_dbl(seq, ~length(.[[1]]))) %>%
  group_by(cluster_grp) %>%
  slice_max(n = 1, order_by = seq_len, with_ties = FALSE) %>%
  ungroup() %>%
  # Optionally filter cluster by exclusion terms
  when(
    !is.null(exclude_terms) ~ filter(., str_detect(cluster, exclude_terms, negate = TRUE)),
    ~ .
  ) %>%
  align_seqs_tbl(name_col = "species")

}


#' Reverse-complement spacer sequences
#'
#' @param plastid_spacers_aligned_trimmed_clusters Trimmed aligned spacer sequences by taxonomic cluster
#'  in a tibble, with alignments in list-column "align_trimmed"
#' @param plastid_spacers_rep_align Status of which spacers need to be reverse-complemented,
#' with one representative sequence per taxonomic cluster. Cluster called "none" are singletons,
#' which do not belong to any clustered alignment.
#'
#' @return Tibble with reverse-complemented sequences. Column "align-trimmed" includes
#' alignments of clusters; column "seq" includes sequences of individual singletons.
#'
reverse_spacers <- function(plastid_spacers_aligned_trimmed_clusters, plastid_spacers_rep_align) {
  
  # Check that input names match arguments
  check_args(match.call())
  
  # Make a tibble showing which clusters were reversed
  clusters_rev_tbl <-
    plastid_spacers_rep_align %>%
    select(cluster, target, reversed) %>%
    filter(cluster != "none") %>%
    unique()
  
  # Reverse-complement clusters
  # will drop any clusters that were excluded when making `clusters_rev_tbl`
  plastid_spacers_aligned_trimmed_rev <-
    plastid_spacers_aligned_trimmed_clusters %>%
    filter(cluster != "none") %>%
    inner_join(clusters_rev_tbl, by = c("target", "cluster")) %>%
    assert(not_na, reversed) %>%
    assert(is.logical, reversed) %>%
    mutate(
      align_trimmed = case_when(
        reversed == TRUE ~ map(align_trimmed, ape::complement),
        reversed == FALSE ~ align_trimmed
      )
    )

  # Extract singletons, but don't rev-complement (MAFFT already did that for us)
  singletons_rev <-
    plastid_spacers_rep_align %>%
    filter(cluster == "none") %>%
    assert(not_na, reversed) %>%
    assert(is.logical, reversed)
  
  # Combine singletons and clusters
  bind_rows(plastid_spacers_aligned_trimmed_rev, singletons_rev) %>%
    assert(not_na, target, cluster, reversed) %>%
    select(-reversed)
  
}

#' Merge subalignments with MAFFT
#' 
#' Modified from ips::mafft.merge(), but allows adding 
# "singleton" seqs that are not part of any subalignment
#'
#' @param subMSA List of sub-alignments to merge; each one a matrix of class "DNAbin"
#' @param other_seqs List of class "DNAbin"; "singleton" sequences that don't belong
#' to any sub-alignment, but should be merged with the sub-alignments.
#' @param method Name of method for MAFFT to use
#' @param gt List of class "phylo"; guide tree
#' @param thread Number of threads to use
#' @param exec Path to MAFFT executable
#' @param quiet Logical; should MAFFT output be printed to screen?
#' @param adjustdirection Logical; should MAFFT attempt to automatically adjust sequence direction?
#'
#' @return Matrix of class "DNAbin"; the merged alignment
#'
mafft_merge <- function (subMSA, other_seqs, method = "auto", gt, thread = -1, exec, quiet = TRUE, adjustdirection = FALSE) 
{
    quiet <- ifelse(quiet, "--quiet", "")
    adjustdirection <- ifelse(adjustdirection, "--adjustdirection", "")
    method <- match.arg(method, c("auto", "localpair", "globalpair", 
        "genafpair", "parttree", "retree 1", "retree 2"))
    method <- paste("--", method, sep = "")
    if (missing(gt)) {
        gt <- ""
    }
    else {
        phylo2mafft(gt)
        gt <- "--treein tree.mafft"
    }
    n <- sapply(subMSA, nrow)
    subMSAtable <- vector(length = length(n))
    init <- 0
    for (i in seq_along(n)) {
        nn <- 1:n[i] + init
        init <- max(nn)
        subMSAtable[i] <- paste(nn, collapse = " ")
    }
    subMSA <- lapply(subMSA, as.list)
    subMSA <- c(subMSA, list(other_seqs))
    subMSA <- do.call(c, subMSA)
    names(subMSA) <- gsub("^.+[.]", "", names(subMSA))
    fns <- vector(length = 3)
    for (i in seq_along(fns)) fns[i] <- tempfile(pattern = "mafft",
        tmpdir = tempdir())
    write(subMSAtable, fns[1])
    ips::write.fas(subMSA, fns[2])
    call.mafft <- paste(exec, method, "--merge", fns[1], quiet,
      adjustdirection,
      gt, "--thread", thread, fns[2], ">", fns[3])
    system(call.mafft, intern = FALSE, ignore.stdout = FALSE)
    res <- length(scan(fns[3], what = "c", quiet = TRUE))
    if (res != 0) {
        res <- ape::read.FASTA(fns[3])
        if (length(unique(sapply(res, length))) == 1) {
            res <- as.matrix(res)
        }
    }
    unlink(fns[file.exists(fns)])
    return(res)
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
#' @param plastome_ncbi_names_raw Species names in plastome data extracted from
#' NCBI taxonomic database
#' @param plastome_metadata Dataframe with column "species" containing
#' original species names from GenBank
#' @param plastome_outgroups Dataframe with metadata on outgroups
#' @param ref_names_parsed Dataframe with parsed scientific names for taxonomic name resolution
#' of fern species; output of ts_parse_names()
#' @param ref_names_data Reference data for taxonomic name resolution of fern species
#' extracted from Catalog of Life data; output of extract_fow_from_col()
#' 
#' @return Dataframe; plastome_metadata with new column `accepted_name` and `species`
#' containing the standardized name attached to it, also a column `outgroup` indicating
#' if the data correspond to outgroup or not
#' 
resolve_pterido_plastome_names <- function(plastome_ncbi_names_raw, 
  plastome_metadata_raw, plastome_outgroups, ref_names_parsed, ref_names_data) {
  
  ### outgroups ###
  # Fetch full scientific names for plastome outgroups
  plastome_outgroup_sci_names <-
    plastome_metadata_raw %>%
    inner_join(select(plastome_outgroups, accession), by = "accession") %>%
    left_join(plastome_ncbi_names_raw, by = "taxid") %>%
    filter(accepted == TRUE) %>%
    select(-accepted) %>%
    # Add one scientific name currently missing
    mutate(
      scientific_name = case_when(
        species == "Sciadopitys verticillata" ~ "Sciadopitys verticillata (Thunb.) Siebold & Zucc.",
        TRUE ~ scientific_name
      )
    ) %>%
    assert(not_na, scientific_name)

  ### ferns ###
  plastome_names_query <- 
    # Exclude outgroups
    plastome_metadata_raw %>%
    anti_join(
      unique(select(plastome_outgroups, accession)),
      by = "accession") %>%
    left_join(plastome_ncbi_names_raw, by = "taxid") %>%
    # Use only accepted name
    filter(accepted == TRUE) %>%
    # Clean names: removes brackets and years
    clean_ncbi_names() %>%
    # (needed `accepted` during cleaning, now ok to drop)
    select(-accepted) %>%
    # Remove any names not identified to species
    filter(str_detect(species, " sp\\. ", negate = TRUE)) %>%
    # Fix some names NCBI got wrong
    mutate(
      scientific_name = case_when(
        # NCBI mistakenly used Drynaria parishii (Bedd.) C.Chr. & Tardieu with basionym as "Pleopeltis parishii Bedd."
        # but the real basionym for Drynaria parishii (Bedd.) C.Chr. & Tardieu is Meniscium parishii Bedd. -> Grypothrix parishii #nolint
        # Pleopeltis parishii Bedd. should point to Drynaria parishii (Bedd.) Bedd. #nolint
        taxid == "2836282" ~ "Drynaria parishii (Bedd.) Bedd.", 
        TRUE ~ scientific_name
      )
    ) %>%
    # Specify query name: scientific name if available, species if not
    mutate(query_name = coalesce(scientific_name, species)) %>%
    assert(not_na, query_name) %>%
    assert(is_uniq, accession)

  # Match names to world ferns
  match_results_plastome <- ts_match_names(
    query = unique(plastome_names_query$query_name),
    reference = ref_names_parsed,
    max_dist = 5, match_no_auth = TRUE, 
    match_canon = TRUE, collapse_infra = TRUE)

  # Resolve synonyms
  match_results_plastome_resolved <- ts_resolve_names(
    match_results_plastome, ref_names_data) %>%
    assert(not_na, everything()) %>%
    # Make sure no names detected by fuzzy match
    # (all names should have already been inspected)
    verify(
      !any(str_detect(match_type, "fuzzy|no_match")),
      error_fun = err_msg(
        "Fuzzy or no match names detected during plastome name resolution")
    )

  ### Combine results ###
  plastome_metadata_ferns_resolved <-
    plastome_metadata_raw %>%
    # filter to only ingroup
    anti_join(
      unique(select(plastome_outgroups, accession)), by = "accession") %>%
    # add queried taxonomic name
    left_join(
      unique(select(plastome_names_query, taxid, query_name)), by = "taxid") %>%
    # add resolved name
    left_join(
      unique(select(
        match_results_plastome_resolved, 
        query_name = query, resolved_name)),
      by = "query_name") %>%
    select(-query_name) %>%
    mutate(outgroup = FALSE) %>%
    assert(is_uniq, accession)

  plastome_metadata_outgroups_resolved <-
    plastome_metadata_raw %>%
    inner_join(
      unique(select(plastome_outgroups, accession)), by = "accession") %>% 
    select(-matches("species|variety")) %>%
    left_join(
      select(
        plastome_outgroup_sci_names, taxid,
        resolved_name = scientific_name),
      by = "taxid") %>%
    mutate(outgroup = TRUE)

  # Combine results
  bind_rows(
    plastome_metadata_ferns_resolved,
    plastome_metadata_outgroups_resolved
  ) %>% 
  # Make sure all accessions are included
  verify(all(accession %in% plastome_metadata_raw$accession)) %>%
  verify(all(plastome_metadata_raw$accession %in% accession)) %>%
  # Drop any without resolved species (only identified to sp., etc)
  filter(!is.na(resolved_name)) %>%
  # Parse resolved name to just species (drops infrasp taxon)
  mutate(
    gn_parse_tidy_quiet(resolved_name) %>% 
      select(taxon = canonicalsimple)
  ) %>%
  mutate(taxon = str_replace_all(taxon, " ", "_")) %>%
  separate(taxon,
    into = c("genus", "sp_epithet", "infrasp_epithet"),
    sep = "_", remove = FALSE, fill = "right") %>%
    assert(not_na, genus, sp_epithet) %>%
    mutate(species = paste(genus, sp_epithet, sep = "_")) %>%
  select(
    sci_name = resolved_name, species, accession,
    subtype, subname, slen, outgroup)
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

#' Merge spacer sub alignments with mafft
#'
#' @param plastid_spacers_aligned_trimmed Tibble with one row per subalignment
#' for clusters, and one row per sequence for singletons
#' @param n_threads Number of threads for mafft to use
#'
#' @return Tibble with list-column `align_trimmed` including alignment;
#' column `cluster` will have the value "combined"
#' 
merge_spacer_alignments <- function(plastid_spacers_reversed, target_select, n_threads) {
  # Extract sub alignments
  sub_alns <- plastid_spacers_reversed %>%
    filter(target == target_select) %>%
    assert(not_na, cluster) %>%
    filter(cluster != "none") %>%
    pull(align_trimmed)
  
  # Extract singletons
  singletons <- plastid_spacers_reversed %>%
    filter(target == target_select) %>%
    assert(not_na, cluster) %>%
    filter(cluster == "none") %>%
    seqtbl_to_dnabin(name_col = "species")
  
  # Merge sub alignments and singletons
  alignment <-
    mafft_merge(
    subMSA = sub_alns,
    other_seqs = singletons,
      thread = n_threads,
      quiet = TRUE,
      method = "retree 2",
      adjustdirection = FALSE,
      exec = "/usr/bin/mafft")
  
  plastid_spacers_reversed %>%
    filter(target == target_select) %>%
    select(target) %>%
    unique() %>%
    mutate(
      cluster = "combined",
      align_trimmed = list(alignment)
    )
}

#' Concatenate plastid and sanger genes
#'
#' For plastome dataset, will include only species with whole plastome data
#' 
#' For Sanger dataset, will include on the genes in `target_loci`
#'
#' @param plastid_genes_aligned_trimmed Tibble with list-column of trimmed
#' DNA alignments; each row is a target gene.
#' @param plastid_spacers_aligned_trimmed Tibble with list-column of trimmed
#' DNA alignments; each row is a target spacer.
#' @param target_loci Character vector of target Sanger loci.
#' @param type_select Type of dataset to filter to: "sanger" or "plastome".
#'
#' @return Tibble with list-column of trimmed DNA alignments (either Sanger
#' or plastome data)
#'
concatenate_plastid_sanger <- function(
  plastid_genes_aligned_trimmed,
  plastid_spacers_aligned_trimmed, target_loci,
  type_select) {
  
  # Make list of plastome species:
  # those in genes that are not target Sanger loci
  plastome_species <-
    plastid_genes_aligned_trimmed %>%
    filter(!target %in% target_loci) %>%
    mutate(species = map(align_trimmed, rownames)) %>%
    select(species) %>%
    unnest(species) %>%
    unique() %>%
    pull(species)
  
  # Filter to plastome alignments by species
  if (type_select == "plastome") {
    res <- bind_rows(
      plastid_genes_aligned_trimmed,
      plastid_spacers_aligned_trimmed) %>%
      select(-cluster) %>%
      mutate(
        align_trimmed = map(
          align_trimmed,
          ~magrittr::extract(., rownames(.) %in% plastome_species, ))
      )
    # Filter to Sanger alignments by gene name
  } else if (type_select == "sanger") {
    res <- bind_rows(
      plastid_genes_aligned_trimmed,
      plastid_spacers_aligned_trimmed) %>%
      select(-cluster) %>%
      filter(target %in% target_loci)
  } else (stop("Must choose 'plastome' or 'sanger' for 'type_select'"))
  
  res
}

#' Concatenate sequences into ape format
#'
#' @param aln_tbl Tibble with list-column of DNA alignments.
#' @param aln_col Name of column with DNA alignments.
#'
#' @return List of class "DNAbin"
#'
concatenate_to_ape <- function(aln_tbl, aln_col = "align_trimmed") {
  do.call(
    ape::cbind.DNAbin,
    c(aln_tbl[[aln_col]],
    fill.with.gaps = TRUE)
  )
}

# Check gene trees ----

#' Build gene trees
#'
#' @param gene_tree_alignment_df 1-row tibble with alignment in list-column
#' `align_trimmed`
#' @param program Name of program used to build tree; "iqtree" or 'fasttree'
#'
#' @return Tibble with tree in list-column `tree`
#' 
build_tree_from_alignment_df <- function(gene_tree_alignment_df, program = "fasttree") {

  assertthat::assert_that(
    program %in% c("fasttree", "iqtree"),
    msg = "Must choose 'fastree' or 'iqtree' for 'program'")
	
	# Extract alignment from df (df should only be one row)
	alignment <- 
		gene_tree_alignment_df %>%
		pull(align_trimmed) %>%
		magrittr::extract2(1)
	
	# Create temp dir for writing the tree (only needed for iqtree)
	temp_dir <- fs::path(tempdir(), digest::digest(alignment))
	
	if(fs::dir_exists(temp_dir)) fs::dir_delete(temp_dir)
	
	fs::dir_create(temp_dir)
	
	# Run single iqtree analysis (no bootstrap)
	if (program == "iqtree") { 
    tree <- iqtree(
		  alignment, wd = temp_dir, 
		  m = "GTR+I+G", nt = "AUTO", 
		  echo = FALSE, 
		  tree_path = fs::path(temp_dir, "alignment.phy.treefile")) 
    } else {
    tree <- fasttree(
		  seqs = alignment,
      mol_type = "dna",
      model = "gtr",
      gamma = FALSE,
		  echo = FALSE, 
		  tree_path = fs::path(temp_dir, "alignment.phy.treefile")) 
  }
	# Cleanup
	if(fs::dir_exists(temp_dir)) fs::dir_delete(temp_dir)
	
	# Return tree in data frame
	gene_tree_alignment_df %>%
		select(-align_trimmed) %>%
		mutate(tree = list(tree))
}

# Dating with treePL ----

#' Get the best smoothing parameter for treePL
#'
#' @param cv_results Results of running treePL with cross-validation to
#' determine optimal rate-smoothing parameter. Output of run_treepl_cv().
#'
#' @return Character
#' 
get_best_smooth <- function(cv_results) {
  # Select the optimum smoothing value (smallest chisq) from cross-validation
  # The raw output looks like this:
  # chisq: (1000) 6.7037e+30
  # chisq: (100) 3673.45
  # etc.
  tibble(cv_result = cv_results) %>%
    mutate(
      smooth = str_match(cv_result, "\\((.*)\\)") %>%
      magrittr::extract(, 2) %>%
      parse_number()) %>%
    mutate(
      chisq = str_match(cv_result, "\\) (.*)$") %>%
      magrittr::extract(, 2) %>%
      parse_number()) %>%
    arrange(chisq) %>%
    slice(1) %>%
    pull(smooth)
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
#'
run_treepl_cv <- function (
  phy, alignment, calibration_dates, 
  write_tree = TRUE,
  cvstart = "1000", cvstop = "0.1",
  cvsimaniter = "5000", 
  plsimaniter = "5000",
  nthreads = "1",
  seed,
  thorough = TRUE, wd) {
  
  # Check that all taxa are in tree
  taxa <- c(calibration_dates$taxon_1, calibration_dates$taxon_2) %>%
    unique()
  
  assertthat::assert_that(all(taxa %in% phy$tip.label),
                          msg = glue(
                            "Taxa in calibration dates not present in tree: 
                            {taxa[!taxa %in% phy$tip.label]}"))
  
  # Write tree to wd
  phy_path <- "undated.tre"
  ape::write.tree(phy, fs::path(wd, phy_path))
  
  # Get number of sites in alignment
  num_sites <- dim(alignment)[2] # nolint
  
  outfile_path <- "treepl_cv_out.txt" # nolint
  
  # Write config file to working directory
  treepl_config <- c(
    glue("treefile = {phy_path}"),
    glue("numsites = {num_sites}"),
    calibration_dates$mrca,
    calibration_dates %>% filter(!is.na(min)) %>% pull(min),
    calibration_dates %>% filter(!is.na(max)) %>% pull(max),
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
  
  config_file_name <- "treepl_cv_config.txt"
  
  readr::write_lines(treepl_config, fs::path(wd, config_file_name))
  
  # Run treePL
  processx::run(
    "treePL", config_file_name, wd = wd,
    stdout = fs::path(wd, "treepl_cv.stdout"),
    stderr = fs::path(wd, "treepl_cv.stderr")
  )
  
  # Return cross-validation results
  readr::read_lines(fs::path(wd, outfile_path))
  
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
  
  # Write tree to wd
  phy_path <- "undated.tre"
  ape::write.tree(phy, fs::path(wd, phy_path))
  
  # Get number of sites in alignment
  num_sites <- dim(alignment)[2] #nolint
  
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
    calibration_dates %>% filter(!is.na(min)) %>% pull(min),
    calibration_dates %>% filter(!is.na(max)) %>% pull(max),
    glue("cvsimaniter = {cvsimaniter}"),
    glue("plsimaniter = {plsimaniter}"),
    glue("seed = {seed}"),
    glue("smooth = {best_smooth}"),
    glue("nthreads = {nthreads}"),
    "prime"
  )
  
  if(thorough) treepl_config <- c(treepl_config, "thorough")
  
  config_file_name <- "treepl_prime_config.txt"
  
  readr::write_lines(treepl_config, fs::path(wd, config_file_name))
  
  # Run treePL
  results <- processx::run(
    "treePL", config_file_name, wd = wd,
    stdout = fs::path(wd, "treepl_prime.stdout"),
    stderr = fs::path(wd, "treepl_prime.stderr")
  )
  
  # Return stdout
  readr::read_lines(fs::path(wd, "treepl_prime.stdout"))
  
}

#' Run treePL
#' 
#' Requires treepl to be installed and on PATH.
#'
#' For more details on treepl options, see
#' https://github.com/blackrim/treePL/wiki
#'
#' @param phy List of class "phylo"; phylogeny
#' @param treepl_config Character vector: treePL config. If provided,
#' will be used in favor of `alignment`, `calibration_dates`, `priming_results`,
#' and `cv_results` to run treePL. 
#' The name of the infile must be "undated.tre".
#' The name of the outfile must be "dated.tre".
#' @param alignment List of class "DNAbin"; alignment
#' @param calibration_dates Calibration points read in with
#' load_calibration_dates()
#' @param priming_results Results of running treePL with `prime` option to
#' determine optional parameters. Output of run_treepl_prime().
#' @param cv_results Results of running treePL with cross-validation to
#' determine optimal rate-smoothing parameter. Output of run_treepl_cv().
#' @param cvsimaniter Start the number of cross validation simulated annealing
#' iterations, default = 5000 for cross-validation
#' @param plsimaniter the number of penalized likelihood simulated annealing
#' iterations, default = 5000
#' @param nthreads Number of threads for treePL to use
#' @param seed Seed for random number generator
#' @param thorough Logical; should the "thorough" setting in
#' treePL be used?
#' @param wd Working directory to run all treepl analyses. If not provided,
#' a temporary one will be created automatically and deleted at the end of the
#' run.
#' The input tree will be written here as "undated.tre".
#' The config fill will be written here as "treepl_config.txt".
#' @param echo Logical; should the output be printed to the screen?
#'
run_treepl <- function(
  phy,
  treepl_config = NULL,
  alignment = NULL,
  calibration_dates = NULL,
  priming_results = NULL,
  cv_results = NULL,
  cvsimaniter = 5000,
  plsimaniter = 5000,
  nthreads = 1,
  seed = 1,
  wd = NULL,
  thorough = TRUE, echo = FALSE) {

  # Save original arg input to wd for checking on this later
  wd_arg <- wd

  # Set up temporary working directory: unique WD for each seed
  if (is.null(wd_arg)) {
    wd <- fs::path(tempdir(), glue::glue("treepl_{seed}"))
    assertthat::assert_that(
      !fs::dir_exists(wd),
      msg = "Temporary treepl directory already exists")
    fs::dir_create(wd)
  }

  # Check that all taxa are in tree
  if(is.null(treepl_config)) {
    taxa <- c(calibration_dates$taxon_1, calibration_dates$taxon_2) %>%
      unique()
    assertthat::assert_that(all(taxa %in% phy$tip.label),
      msg = glue(
        "Taxa in calibration dates not present in tree:
        {taxa[!taxa %in% phy$tip.label]}"))
  }

  # Write tree to wd
  phy_path <- "undated.tre"
  ape::write.tree(phy, fs::path(wd, phy_path))

  # If the treePL config is provided, run that instead
  if (!is.null(treepl_config)) {
    readr::write_lines(treepl_config, fs::path(wd, "treepl_config.txt"))
    # Run treePL
    processx::run("treePL", "treepl_config.txt", wd = wd, echo = echo)
    res <- ape::read.tree(fs::path(wd, "dated.tre"))
    return(res)
  }

  # Get number of sites in alignment
  num_sites <- dim(alignment)[2] #nolint

  # Get best smoothing parameter
  best_smooth <- get_best_smooth(cv_results)

  # Set name of output file
  outfile_path <- "dated.tre"

  # Write config file to working directory
  treepl_config <- c(
    glue("treefile = {phy_path}"),
    glue("numsites = {num_sites}"),
    glue("smooth = {best_smooth}"),
    calibration_dates$mrca,
    calibration_dates %>% filter(!is.na(min)) %>% pull(min),
    calibration_dates %>% filter(!is.na(max)) %>% pull(max),
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

  if (thorough) treepl_config <- c(treepl_config, "thorough")

  config_file_name <- "treepl_config.txt"

  readr::write_lines(treepl_config, fs::path(wd, config_file_name))

  # Run treePL
  processx::run("treePL", config_file_name, wd = wd, echo = echo)

  # Read in tree
  dated_tree <- ape::read.tree(fs::path(wd, outfile_path))

  # Delete any temporary wd
  if (fs::dir_exists(wd) && is.null(wd_arg)) fs::dir_delete(wd)

  return(dated_tree)

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
  
  label <- paste3(as.character(length), units, sep = " ")
  
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

  require(ape) # for combining seqs with c()
  
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
    do.call(c, .)
  
  # Write out pseudo reference genome
  ape::write.FASTA(pseudo_ref_genome, fs::path(temp_dir, "ref.fasta"))
  
  # Load all the sorted reads from HybPiper for this sample, combine
  short_reads <- list.files(
    paste0("intermediates/hybpiper/", sample), 
    pattern = "interleaved.fasta", 
    full.names = TRUE, recursive = TRUE) %>%
    map(ape::read.FASTA) %>%
    do.call(c, .)
  
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

#' Make accessions table in long format
#' 
#' (one row per accession)
#'
#' @param raw_meta Tibble; Raw metadata for GenBank accessions in Sanger dataset
#' @param sanger_seqs_combined_filtered Tibble; DNA sequences and metadata for
#' Sanger dataset after removing rogues by all-by-all BLAST
#' @param plastome_seqs_combined_filtered Tibble; DNA sequences and metadata for
#' plastome dataset, filtered to single best accession per species
#' @param ncbi_names_query Tibble; NCBI species names to query against 
#' taxonomic standard
#' @param sanger_accessions_selection Tibble; Final selection of accessions
#' for Sanger dataset
#' @param plastome_metadata_renamed Tibble; Metadata for plastome dataset
#' renamed after resolving species names against taxonomic standard
#' @param plastome_metadata_raw Tibble; Raw metadata for GenBank accessions 
#' in plastome dataset
#' @param plastome_ncbi_names_raw Tibble; Species names in plastome data 
#' extracted from NCBI taxonomy
#'
#' @return Tibble with one row per accession including species (tips in FTOL),
#' scientific name, accession length, etc.
#' 
make_long_acc_table <- function(
  raw_meta, sanger_seqs_combined_filtered,
  plastome_seqs_combined_filtered,
  ncbi_names_query, sanger_accessions_selection,
  plastome_metadata_renamed,
  plastome_metadata_raw,
  plastome_ncbi_names_raw) {

  # Check argument names
  check_args(match.call())
  
  ## Plastome data ##
  # Tibble of NCBI accepted names
  plastome_ncbi_accepted_names <-
    plastome_ncbi_names_raw %>%
    filter(accepted == TRUE) %>%
    mutate(
      ncbi_name = coalesce(scientific_name, species) %>%
        # Drop years
        str_remove_all(", [0-9]+$")) %>%
    select(taxid, ncbi_name) %>%
    assert(not_na, everything()) %>%
    assert(is_uniq, everything())
  
  # Tibble of taxid + accession
  plastome_taxid_acc <-
    plastome_metadata_raw %>%
    select(taxid, accession) %>%
    unique() %>%
    assert(not_na, everything()) 
  
  # Tibble of plastome seqlengths
  # (actual seq length of non-missing bases in final dataset)
  plastome_acc_seqlen <-
    plastome_seqs_combined_filtered %>%
    select(accession, target, seq) %>%
    mutate(seq_len = map_dbl(seq, ~count_non_missing(.[1]))) %>%
    group_by(accession) %>%
    summarize(seq_len = sum(seq_len))
  
  # Tibble of plastome data in long format
  plastome_data <-
    plastome_seqs_combined_filtered %>%
    select(species, accession) %>%
    unique() %>%
    assert(not_na, everything()) %>%
    assert(is_uniq, everything()) %>%
    mutate(locus = "plastome") %>%
    left_join(plastome_metadata_renamed, by = c("species", "accession")) %>%
    select(species, accession, locus, sci_name, outgroup) %>%
    # Add seq len
    left_join(plastome_acc_seqlen, by = "accession") %>%
    # Add taxid
    left_join(plastome_taxid_acc, by = "accession") %>%
    # Add NCBI name
    left_join(plastome_ncbi_accepted_names, by = "taxid") %>%
    assert(not_na, everything()) %>%
    assert(is_uniq, species, accession)
  
  ## Sanger data
  
  # Tibble of taxid + scientific name
  # unique taxid, but not resolved_name
  taxid_sci_name <-
    sanger_seqs_combined_filtered %>%
    select(taxid, sci_name = resolved_name) %>%
    unique() %>%
    assert(not_na, everything()) %>%
    assert(is_uniq, taxid) 
  
  # Tibble of taxid + accepted NCBI name
  ncbi_accepted_names <-
    ncbi_names_query %>%
    filter(accepted == TRUE) %>%
    mutate(ncbi_name = coalesce(scientific_name, species)) %>%
    select(taxid, ncbi_name) %>%
    assert(not_na, everything()) %>%
    assert(is_uniq, everything())
  
  # Tibble of taxid, accession, and locus (gene)
  # combination of accession + locus is unique
  taxid_accession_locus <-
    raw_meta %>%
    select(taxid, accession, locus = target) %>%
    unique() %>%
    assert(not_na, everything()) %>%
    assert_rows(col_concat, is_uniq, locus, accession) 
  
  # Convert Sanger seqs to long format
  sanger_seq_len <-
    sanger_accessions_selection %>%
    select(species, matches("seq_len")) %>%
    select(-total_seq_len) %>%
    pivot_longer(names_to = "locus", values_to = "seq_len", -species) %>%
    mutate(locus = str_remove_all(locus, "seq_len_")) %>%
    filter(!is.na(seq_len)) %>%
    filter(seq_len > 0)
  
  ## Combine data into final table
  # Convert Sanger seqs to long format
  sanger_accessions_selection %>%
    select(species, matches("accession")) %>%
    pivot_longer(names_to = "locus", values_to = "accession", -species) %>%
    mutate(locus = str_remove_all(locus, "accession_")) %>%
    filter(!is.na(accession)) %>%
    # Add taxid by accession + locus
    left_join(taxid_accession_locus, by = c("accession", "locus")) %>%
    # Add sci name by taxid
    left_join(taxid_sci_name, by = "taxid") %>%
    left_join(ncbi_accepted_names, by = "taxid") %>%
    # Add seq len
    left_join(sanger_seq_len, by = c("species", "locus")) %>%
    assert(not_na, everything()) %>%
    assert_rows(col_concat, is_uniq, locus, accession) %>%
    mutate(outgroup = FALSE) %>%
    # Remove species in plastome data
    anti_join(plastome_data, by = "species") %>%
    # Add plastome species
    bind_rows(plastome_data) %>%
    arrange(outgroup, species, locus) %>%
    select(
      species, locus, accession, seq_len, 
      sci_name, ncbi_name, ncbi_taxid = taxid,
      outgroup) %>%
    assert(not_na, everything()) %>%
    assert_rows(col_concat, is_uniq, species, locus)
}

#' Make accessions table in wide format
#'
#' (one row per species)
#'
#' @param acc_table_long Tibble; accession data in long format.
#' @param sanger_accessions_selection Metadata for Sanger fern accessions.
#'
#' @return Tibble; accession data in wide format.
#'
make_wide_acc_table <- function(acc_table_long, sanger_accessions_selection) {
  check_args(match.call())

  # Make tibble of plastome species
  species_plastome <-
    acc_table_long %>%
    filter(locus == "plastome") %>%
    unique()

  # Read in data on join method, voucher, and publication (Sanger ferns only)
  species_join_voucher_pub <-
    sanger_accessions_selection %>%
    select(species, join_by, specimen_voucher, publication) %>%
    anti_join(species_plastome, by = "species")

  # Make tibble of species + outgroups
  species_og <-
    acc_table_long %>%
    select(species, outgroup) %>%
    unique()

  # Assemble accessions in wide format
  acc_table_long %>%
    select(species, locus, accession) %>%
    tidyr::pivot_wider(names_from = locus, values_from = accession) %>%
    select(
      species, atpA, atpB, matK, rbcL, rps4,
      `rps4-trnS`, `trnL-trnF`, plastome) %>%
    left_join(species_join_voucher_pub, by = "species") %>%
    left_join(species_og, by = "species") %>%
    assert(is_uniq, species) %>%
    assert(not_na, species) %>%
    arrange(outgroup, species)
}

# Format data for ftolr ----

#' Make a partition (by locus) table for a set of concatenated
#' DNA alignments
#'
#' @param aln_tbl Tibble with list-column of DNA alignments.
#' @param aln_seq Matrix of class "DNAbin": the concatenated DNA alignments.
#' (only used for double-checking that aln_tbl and aln_seq are have identical
#' sequences)
#'
#' @return Tibble with columns "locus", "start", "end"
#'
make_parts_table <- function(aln_tbl, aln_seq) {
  # Get start and end position of each gene in concatenated alignment
  res <- aln_tbl %>%
  mutate(
    nbp = map_dbl(align_trimmed, ncol),
    end = cumsum(nbp),
    start = end - nbp + 1) %>%
  select(locus = target, start, end) %>%
  assert(is_uniq, everything()) %>%
  assert(not_na, everything())

  # Double check that aln_tbl and aln_seq (alignment actually used
  # for phy analysis) match
  assertthat::assert_that(
    isTRUE(all.equal(
      aln_seq,
      concatenate_to_ape(aln_tbl)
    )),
  msg = "aln_tbl and aln_seq don't match"
  )

  return(res)
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

# Names to exclude when reading in NCBI taxonomic database
ncbi_db_names_to_exclude <- function() {
  c(
    "Archangiopteris hokouensis Ching, 1958, non Angiopteris hokouensis Ching, 1959", #nolint
    "Polystichum imbricans subsp. curtum (Ewan) D.H.Wagner, 1979")
}

plastome_ncbi_db_names_to_exclude <- function() {
  # Superfluous with Danaea sellowiana C.Presl, 1845
  "Danaea sellowiana Pr.in Corda., 1845"
}

#' Extract taxonomic names from an NCBI taxonomy database dump file
#'
#' Taxonomy dump files can be downloaded from the NCBI FTP server
#' https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/
#' The original data format should be used (not "new")
#'
#' There is only one accepted name per taxid
#'
#' @param taxdump_zip_file Path to zip file with NCBI taxonomy database; must
#' contain a file called "names.dmp".
#' @param taxid_keep Dataframe (tibble) with a column 'taxid' including the NCBI
#' taxids of the names to be extracted.
#' @param names_exclude Character vector; names to exclude from results
#' (optinal).
#' @param workers Number of cores to use in parallel during processing.
#'
#' @return Tibble with columns with columns "taxid" (NCBI taxonomy database ID),
#' "species" (species name without author), "accepted" (logical indicated if
#' name is accepted name or not), "scientific_name" (species name with author)
#'
extract_ncbi_names <- function(taxdump_zip_file, taxid_keep, names_exclude = NULL, workers = 1) {
   # Unzip names.dmp to a temporary directory
   temp_dir <- tempdir(check = TRUE)

   utils::unzip(
      taxdump_zip_file, files = "names.dmp",
      overwrite = TRUE, junkpaths = TRUE, exdir = temp_dir)

   # Load raw NCBI data
   ncbi_raw <-
      fs::path(temp_dir, "names.dmp") %>%
      readr::read_delim(
        delim = "\t|\t", col_names = FALSE,
        col_types = cols(.default = col_character())
      )

   # Delete temporary unzipped file
   fs::file_delete(fs::path(temp_dir, "names.dmp"))

   # Prune raw NBCI names to names in metadata
   ncbi_names <-
      ncbi_raw %>%
      # Select only needed columns
      transmute(
         taxid = as.character(X1),
         name = X2,
         class = X4) %>%
      # Make sure all taxids from metadata are in NCBI data
      verify(all(taxid_keep$taxid %in% taxid)) %>%
      # Filter to only taxids in metadata
      inner_join(
         unique(select(taxid_keep, taxid)),
         by = "taxid"
      ) %>%
      # Make sure there are no hidden fields in `class`
      verify(all(str_count(class, "\\|") == 1)) %>%
      # Drop field separators in `class`
      mutate(class = str_remove_all(class, "\\\t\\|")) %>%
      # Only keep useful names: exclude common names,
      # alternative spellings (`equivalent name`), type material,
      # temporary names with 'sp.' (`includes`)
      filter(class %in% c("authority", "scientific name", "synonym"))

   # Optionally exclude some names
   if(!is.null(names_exclude)) {
      ncbi_names <-
         ncbi_names %>%
         # Make sure all taxids from metadata are in NCBI data
         verify(all(names_exclude %in% .$name)) %>%
         # Exclude some superfluous sci names: these cause multiple accepted names 
         # for a given taxid
      filter(!name %in% names_exclude)
   }

   # Parse names
   # A little slow, so do in parallel.
   # Set backend for parallelization
   if (workers > 1) future::plan(future::multisession, workers = workers)

   ncbi_names_parsed <-
      ncbi_names %>%
      group_by(taxid) %>%
      nest() %>%
      ungroup() %>%
      mutate(
         data = map(data, parse_ncbi_tax_record)
      ) %>%
      unnest(data)

   # Close parallel workers
   future::plan(future::sequential)

   ncbi_names_parsed %>%
     # Make sure all parsed names have only one accepted name
     group_by(taxid) %>%
     mutate(n_accepted = sum(accepted)) %>%
     ungroup() %>%
     verify(all((n_accepted) == 1)) %>%
     select(-n_accepted)
}

#' Extract the first target word from a string
#'
#' The first occurrence of one of the target words
#' will be extracted from the string, if there is a match.
#' Only exact matching is used (no grep expressions)
#'
#' @param string Character vector of length 1.
#' @param target Character vector.
#'
#' @return Character vector of length 1, or NA if no hits.
#' @examples
#' str_extract_first_target_single("My [name] is", c("[name]", "is"))
str_extract_first_target_single <- function(string, target) { #nolint
   hits <- purrr::map_chr(
      target,
      ~stringr::str_extract(string, stringr::fixed(.)))
   if (all(is.na(hits))) return(NA)
   hits[!is.na(hits)][[1]]
}

#' Extract the first target substring from a string
#'
#' Vectorized version of str_extract_first_target_single()
#'
#' @param string Character vector. String to extract substrings from.
#' @param target Character vector. Target substrings to extract.
#'
#' @return Character vector of length 1, or NA if no hits.
#' @examples
#' str_extract_first_target(
#'   c("My [name] is", "[bar] foo", "a"),
#'   c("[bar]", "is")
#' )
str_extract_first_target <- function(string, target) {
   purrr::map_chr(string, ~str_extract_first_target_single(., target = target))
}

#' Parse a single record from the NCBI taxonomy database
#'
#' The NCBI taxonomy database contains names in the `names.dmp` file.
#' A single record (corresponding to one `taxid`) looks like this:
#'
#' # A tibble: 4 × 3
#'   taxid  name                                               class
#'   <chr>  <chr>                                              <chr>
#' 1 857989 Alansmia glandulifera (A.Rojas) Moguel & M.Kessler authority
#' 2 857989 Alansmia glandulifera                              scientific name
#' 3 857989 Terpsichore glandulifera A.Rojas                   authority
#' 4 857989 Terpsichore glandulifera                           synonym
#'
#' @param record Tibble (dataframe); a single record of NCBI taxonomy data
#'
#' @return Tibble; parsed data with columns "species", "accepted", and
#'   "scientific_name"
#' @examples
#' record <- tribble(
#'   ~taxid, ~name, ~class,
#'   "857989", "Alansmia glandulifera (A.Rojas) Moguel & M.Kessler", "authority", #nolint
#'   "857989", "Alansmia glandulifera", "scientific name",
#'   "857989", "Terpsichore glandulifera A.Rojas", "authority",
#'   "857989", "Terpsichore glandulifera", "synonym"
#' )
#' parse_ncbi_tax_record(record)
parse_ncbi_tax_record <- function(record) {

   # Each record should have exactly 1 accepted species name
   accepted_species <- record %>%
      filter(class == "scientific name") %>%
      pull(name)
   assertthat::assert_that(
      length(accepted_species) == 1,
      msg = "Not exactly 1 accepted scientific name detected")
   # Confusingly, "synonyms" may sometimes contain the accepted name :/
   synonyms <- record %>%
      filter(class == "synonym") %>%
      pull(name)
   # Scientific names (with author) have the class "authority"
   sci_names <- record %>%
      filter(class == "authority") %>%
      pull(name)

   # The results should always have at least taxon ID and species
   species_dat <- tibble(species = accepted_species, accepted = TRUE)

   # Create empty tibbles to hold other name data
   acc_sci_names_dat <- tibble()
   syn_sci_names_dat <- tibble()

   # If other sci names are given, one of them should be the species name
   if (length(sci_names) > 0)
      acc_sci_names_dat <- tibble(scientific_name = sci_names) %>%
      mutate(species = str_extract(
         scientific_name, stringr::fixed(accepted_species))) %>%
      filter(!is.na(species)) %>%
      mutate(accepted = TRUE)

   # If "synonym" and other sci names are given, one (or more) of them are
   # the synonym
   if (length(synonyms) > 0 && length(sci_names) > 0)
      syn_sci_names_dat <- tibble(scientific_name = sci_names) %>%
      # Each sci name should only correspond to max. one synonym
      mutate(species = str_extract_first_target(scientific_name, synonyms)) %>%
      filter(!is.na(species)) %>%
      mutate(accepted = FALSE)

   # If "synonym" is present but no other sci names are given, "synonym" is
   # actually the scientific name of the species
   if (length(synonyms) > 0 && length(sci_names) == 0)
      acc_sci_names_dat <- tibble(scientific_name = synonyms) %>%
      mutate(species = str_extract(
         scientific_name, stringr::fixed(accepted_species))) %>%
      filter(!is.na(species)) %>%
      mutate(accepted = TRUE)

   # Combine scientific names of synonyms and accepted names
   combined_sci_names_dat <- bind_rows(syn_sci_names_dat, acc_sci_names_dat)

   # Join to accepted species with taxon ID
   if (nrow(combined_sci_names_dat) > 0)
      species_dat <- full_join(
         species_dat, combined_sci_names_dat, by = c("species", "accepted"))

   species_dat
}

#' Clean up species names extracted from NCBI taxonomy database
#'
#' @param ncbi_names_raw Tibble with columns `taxid` `species` `accepted` 
#' and `scientific_name`.
#'
#' @return Tibble
#' 
clean_ncbi_names <- function(ncbi_names_raw) {
  ncbi_names_raw %>%
    mutate(
      # Remove brackets around species name
      # (notation in NCBI taxonomic db that genus level taxonomy is uncertain)
      species = str_remove_all(species, "\\[|\\]"),
      # Remove year after authorship in scientific name
      # not all years entered with four digits, so match >1 digit
      scientific_name = str_remove_all(scientific_name, ", [0-9]+"),
      # Fix some scientific names
      scientific_name = case_when(
        # missing author
        species == "Dryopteris basisora" ~ "Dryopteris basisora Christ",
        TRUE ~ scientific_name
      ),
      # Drop "non" part of name,
      # otherwise might match on the "non" author!
      scientific_name = str_remove_all(scientific_name, ", non [^$]+$") #nolint
    ) %>%
    # Mixup in NCBI taxonomy: 160848 should be Hymenophyllum baileyanum Domin
    # at least of 2022-02-02
    filter(!(taxid == "160848" & accepted == FALSE)) %>%
      mutate(
        scientific_name = case_when(
          taxid == "160848" ~ "Hymenophyllum baileyanum Domin",
          TRUE ~ scientific_name
        ),
        species = case_when(
          taxid == "160848" ~ "Hymenophyllum baileyanum",
          TRUE ~ species
        )
      ) %>%
      # Mixup in NCBI taxonomy: 295380 should be Hymenophyllum badium Hook. & Grev.
    filter(!(taxid == "295380" & species == "Trichomanes badium")) %>%
    mutate(
      accepted = case_when(
          taxid == "295380" ~ TRUE,
          TRUE ~ accepted
        )
     ) %>%
     # Mixup in NCBI taxonomy: 449813 should be Arthromeris wallichiana (Spreng.) Ching
     # not Selliguea wallichiana Hook. -> Loxogramme wallichiana (Hook.) M.G.Price
    filter(!(taxid == "449813" & scientific_name == "Selliguea wallichiana Hook.")) %>%
    filter(!(taxid == "449813" & scientific_name == "Polypodium wallichianum Spreng.")) %>%
    mutate(
      accepted = case_when(
          taxid == "449813" ~ TRUE,
          TRUE ~ accepted
        )
     ) %>%
     # Mixup: Acrostichum scandens should be Acrostichum scandens Raddi, not
     # Acrostichum scandens Bory ex Fee for accs GU376696 GU376547
     mutate(
       scientific_name = case_when(
          taxid == "861201" ~ "Acrostichum scandens Raddi",
          TRUE ~ scientific_name
        )
     ) %>%
     # Mixup: Cyathea affinis should be Cyathea affinis Brack. (from Samoa), not
     # Cyathea affinis (G.Forst.) Sw. for acc MT657764 
     mutate(
       scientific_name = case_when(
          taxid == "2853751" ~ "Cyathea affinis Brack.",
          TRUE ~ scientific_name
        )
     ) %>%
     # Trichomanes bimarginatum has two ambiguous synonyms,
     # Trichomanes bimarginatum (Bosch) Bosch -> Didymoglossum bimarginatum (Bosch) Ebihara & K. Iwats. #nolint
     # Trichomanes bimarginatum Bedd. -> Vandenboschia birmanica (Bedd.) Ching
     # in this dataset there is only one sequence (AB257494) and it is
     # Trichomanes bimarginatum (Bosch) Bosch
     mutate(
       scientific_name = case_when(
          taxid == "381227" ~ "Trichomanes bimarginatum (Bosch) Bosch",
          TRUE ~ scientific_name
        )
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
  # Exclude names from consideration that aren't fully identified to species, 
  # environmental samples, or hybrid formulas.
  # Hybrid names *can* be parsed:
  # - "Equisetum x ferrissii" (x before specific epithet)
  # - "x Cystocarpium roskamianum" (x before nothogenus)
  # Hybrid formulas *can't* be parsed:
  # - "Cystopteris alpina x Cystopteris fragilis" (x before another species)
  ncbi_names_exclude <-
    ncbi_names %>%
    filter(
      str_detect(species, " sp\\.| aff\\.| cf\\.| × [A-Z]| x [A-Z]|environmental sample") | #nolint
        str_detect(scientific_name, " sp\\.| aff\\.| cf\\.| × [A-Z]| x [A-Z]|environmental sample") | #nolint
        str_count(species, " ") < 1 | 
        str_count(scientific_name, " ") < 1
      )
  
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

#' Select NCBI names for second round of taxonomic name resolution
#' 
#' Names that are considered synonyms by NCBI and have full scientific name (with author)
#'
#' @param match_results_resolved Dataframe; output of ts_resolve_names() 
#' @param ncbi_names Names downloaded from NCBI taxonomy database.
#'
#' @return Dataframe (tibble)
#' 
select_ncbi_names_round_2 <- function(match_results_resolved_round_1, ncbi_names) {
  
  # Get IDs of all resolved names from round 1
  ncbi_id_resolved <-
    match_results_resolved_round_1 %>%
    left_join(ncbi_names, by = c(query = "scientific_name")) %>%
    filter(!is.na(resolved_name)) %>%
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
#' @param match_results_resolved_round_1 Dataframe; output of ts_resolve_names() 
#' @param match_results_resolved_round_2 Dataframe; output of ts_resolve_names() 
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
    filter(!is.na(resolved_name)) %>%
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
#' @param ... Dataframes; output of ts_resolve_names() 
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
    filter(!is.na(resolved_name)) %>% 
    select(taxid, resolved_name) %>%
    unique() %>% 
    add_count(taxid) %>%
    # Drop any taxid with multiple distinct resolved names
    filter(n == 1) %>%
    assert(is_uniq, taxid) %>%
    select(-n) %>%
    # Add taxon (e.g., 'Foogenus barspecies fooinfraspname')
    mutate(
      rgnparser::gn_parse_tidy(resolved_name) %>% 
        select(taxon = canonicalsimple)
    ) %>%
    mutate(taxon = str_replace_all(taxon, " ", "_")) %>%
    separate(taxon, into = c("genus", "sp_epithet", "infrasp_epithet"), sep = "_", remove = FALSE, fill = "right") %>%
    assert(not_na, genus, sp_epithet) %>%
    mutate(species = paste(genus, sp_epithet, sep = "_"))
}

# Phylogenetic analysis ----

#' Run IQ-TREE
#'
#' For details, see http://www.iqtree.org/doc/
#'
#' @param alignment DNA alignment to use for phylogenetic analysis. Must be
#'   matrix (i.e., aligned sequences) of class DNAbin
#' @param aln_path Path to DNA alignment. Either alignment, or aln_path must be
#'   provided, but not both
#' @param tree_path Optional; path to tree(s) written out by IQ-TREE,
#'   useful if this differs from default alignment name.
#' @param wd Path to working directory. The alignment and IQ-TREE intermediate
#'   files and results will be written here.
#' @param bb Optional; number of ultrafast bootstrap replicates to run.
#' @param nt Optional; number of cores to use. Set to "AUTO" to determine
#'   automatically.
#' @param alrt Optional; number of SH-aLRT tests to run.
#' @param m Optional; specify model. If no model is given, ModelTest will be run
#'   to identify the best model for the data.
#' @param redo Logical; should the analysis be redone from scratch if output
#'   from previous runs is present?
#' @param spp Path to partition file.
#' @param seed Optional; Specify a random number seed to reproduce a previous
#'   run.
#' @param echo Logical; should STDERR be written to the screen?
#' @param other_args Other arguments to pass to IQ tree; must be entered as a
#'   character vector with the name and value of each argument separately. For
#'   example, c("-pers", "0.2", "-nstop", "500").
#' @param ... Additional arguments; not used by this function.
#'
#' @return List; either a single phylogenetic tree (list of class "phylo"),
#' or list containing phylogenetic trees
#'
#' @examples
#' \dontrun{
#' library(ape)
#' temp_dir <- fs::dir_create(tempdir(), "blah")
#' data(woodmouse)
#' # Rapid boot-strap tree with 1000 replicates on best-fitting model
#' tree <- iqtree(woodmouse, temp_dir, bb = 1000, echo = TRUE)
#' plot(tree)
#' # Check the optimum number of cores to use for GTR+I+G model
#' iqtree(
#'   woodmouse,
#'   temp_dir,
#'   m = "GTR+I+G", nt = "AUTO", echo = TRUE, redo = TRUE)
#' fs::dir_delete(temp_dir)
#' }
iqtree <- function(alignment = NULL, wd = getwd(),
                   aln_path = NULL,
                   tree_path = NULL,
                   bb = NULL, nt = NULL, alrt = NULL, m = NULL, redo = FALSE,
                   spp = NULL,
                   seed = NULL,
                   echo = FALSE,
                   other_args = NULL, ...) {
  
  assertthat::assert_that(
    !is.null(alignment) | !is.null(aln_path),
    msg = "Either alignment or aln_path must be provided, but not both")
  
  assertthat::assert_that(
    is.null(alignment) | is.null(aln_path),
    msg = "Either alignment or aln_path must be provided, but not both")
  
  assertthat::assert_that(assertthat::is.dir(wd))
  
  assertthat::assert_that(is.logical(echo))
  
  assertthat::assert_that(is.logical(redo))
  
  if (!is.null(bb))
    assertthat::assert_that(assertthat::is.number(bb))
  
  if (!is.null(alrt))
    assertthat::assert_that(assertthat::is.number(alrt))
  
  if (!is.null(nt))
    assertthat::assert_that(
      assertthat::is.number(nt) | assertthat::is.string(nt))
  
  if (!is.null(m))
    assertthat::assert_that(assertthat::is.string(m))
  
  if (!is.null(spp))
    assertthat::assert_that(assertthat::is.readable(spp))
  
  if (!is.null(seed))
    assertthat::assert_that(assertthat::is.number(seed))
  
  wd <- fs::path_norm(wd)
  
  # check that iqtree is installed and on the PATH
  tryCatch({
    processx::run("iqtree2", "-h", echo = FALSE)
  }, warning = function(w) {
    stop("iqtree not installed and on path")
  }, error = function(e) {
    stop("iqtree not installed and on path")
  }, finally = {
    TRUE
  })
  
  # Write alignment to working directory in phylip format if alignment
  # is provided via R as DNAbin
  if (is.null(aln_path)) {
    assertthat::assert_that(inherits(alignment, "DNAbin"),
                            msg = "alignment must be of class 'DNAbin'")
    assertthat::assert_that(
      is.matrix(alignment),
      msg = "alignment must be a matrix (not a list of unaligned sequences)")
    
    aln_path <- fs::path(wd, deparse(substitute(alignment))) %>%
      fs::path_ext_set("phy")
    
    phangorn::write.phyDat(alignment, aln_path, format = "phylip")
  }
  
  assertthat::assert_that(assertthat::is.readable(aln_path))
  
  # Set up arguments
  iqtree_arguments <- c(
    "-s", fs::path_abs(aln_path),
    if (!is.null(bb)) "-bb",
    bb,
    if (!is.null(alrt)) "-alrt",
    alrt,
    if (!is.null(nt)) "-nt",
    nt,
    if (!is.null(m)) "-m",
    m,
    if (!is.null(seed)) "-seed",
    seed,
    if (!is.null(spp)) "-spp",
    fs::path_abs(spp),
    if (isTRUE(redo)) "-redo",
    other_args
  )
  
  # Run iqtree command
  processx::run(
    "iqtree2",
    iqtree_arguments, wd = wd, echo = echo,
    # Include env variable as workaround for initial parsimony analysis
    # using all cores
    # https://github.com/iqtree/iqtree2/issues/18
    env = c("current", OMP_NUM_THREADS = "1"))
  
  # Read in resulting tree(s)
  # Default: use default treefile if tree_path not provided
  if (is.null(tree_path)) {
    tree_path <- fs::path(wd, deparse(substitute(alignment))) %>%
      fs::path_ext_set(".phy.treefile")
  }

  # Return single tree if only one file in tree_path
  if(length(tree_path) == 1) {
    assertthat::assert_that(assertthat::is.readable(tree_path))
    res <- ape::read.tree(tree_path)
  }

  # Return list of trees if multiple files in tree_path
  if(length(tree_path) > 1) {
    # Set up results list to have same
    # names as tree_path
    res <- vector(length = length(tree_path))
    names(res) <- names(tree_path)
    res <- as.list(res)
    for(i in seq_along(tree_path)) {
      assertthat::assert_that(assertthat::is.readable(tree_path[[i]]))
      # although we check if file is readable
      # beware that ape will crash R if it is not a newick file!
      res[[i]] <- ape::read.tree(tree_path[[i]])
    }
  }

  return(res)
  
}

#' Infer a phylogenetic tree using FastTree
#'
#' Required FastTree to be installed and on user's PATH.
#'
#' @param seqs DNA alignment of class DNAbin
#' @param mol_type Molecule type; either "dna" or "aa" (proteins)
#' @param model Model to use for phylogenetic analysis. Choose "wag" (WAG+CAT),
#'   "lg" (LG+CAT), "gtr" (GTR+CAT), "jc" (Jukes-Cantor + CAT)
#' @param gamma Logical; should branch lengths be rescaled and Gamma20-based
#'   likelihood calculated?
#' @param echo Logical; should STDERR and STDOUT be printed to the screen?
#' @param other_args Other arguments to pass to fasttree; must be entered as a
#'   character vector with the name and value of each argument separately.
#' @param ... Additional arguments; not used by this function.
#'
#' @return List of class "phylo".
#' @references http://www.microbesonline.org/fasttree/
#'
#' @examples
#' \dontrun{
#' library(ape)
#' data(woodmouse)
#' fasttree(woodmouse)
#' }
fasttree <- function (seqs, mol_type = "dna", model = "gtr", gamma = FALSE, other_args = NULL, echo = FALSE, ...) {

  # Make sure input types are correct
  assertthat::assert_that(inherits(seqs, "DNAbin"),
                          msg = "seqs must be of class DNAbin")
  assertthat::assert_that(is.matrix(seqs),
                          msg = "seqs must be in matrix format (aligned)")
  assertthat::assert_that(assertthat::is.string(mol_type))
  assertthat::assert_that(mol_type %in% c("dna", "aa"),
                          msg = "mol_type must be either 'dna' or 'aa'")
  assertthat::assert_that(assertthat::is.string(model))
  assertthat::assert_that(model %in% c("wag", "lg", "gtr", "jc"),
                          msg = "model must be either 'wag', 'lg', 'gtr', or 'jc'")
  assertthat::assert_that(is.logical(gamma))
  assertthat::assert_that(is.logical(echo))

  # Write out alignment to temp file
  temp_wd <- tempdir()
  ape::write.FASTA(seqs, fs::path(temp_wd, "seqs.fasta"))

  # Modify arguments for processx::run()
  mol_type <- if (mol_type == "dna") "-nt" else NULL
  model <- if (model %in% c("wag", "lg", "gtr")) paste0("-", model) else NULL
  gamma <- if (isTRUE(gamma)) "-gamma" else NULL
  alignment_file <- fs::path(temp_wd, "seqs.fasta")

  args <- c(mol_type,
            model,
            gamma,
            alignment_file,
            other_args)

  # Run command
  results <- processx::run(
    "fasttree",
    args, wd = temp_wd, echo = echo)

  # Convert tree to ape format by writing out then reading in
  readr::write_file(results$stdout, fs::path(temp_wd, "tre.fasta"))
  ape::read.tree(fs::path(temp_wd, "tre.fasta"))
}

#' Remove node labels from a phylogenetic tree
#' 
#' Used for generating a constraint tree with no bootstrap support values
#'
#' @param phy Phylogenetic tree; list of class "phylo".
#'
#' @return Phylogenetic tree with no node labels
#' 
remove_node_labels <- function(phy) {
  phy$node.label <- NULL
  phy
}

#' Make a single bootstrap tree with IQTREE
#'
#' For making trees with the same topology but different branchlengths
#' so we can obtain a range of age estimates with treePL
#'
#' Either alignment or aln_path should be provided, but not both
#'
#' @param alignment Matrix; DNA alignment.
#' @param aln_path Path to DNA alignment in phylip format.
#' @param constraint_tree Constraint tree; list of class "DNA bin".
#' @param m Optional; specify model. If no model is given, ModelTest will be run
#'   to identify the best model for the data.
#' @param nt Optional; number of cores to use. Set to "AUTO" to determine
#'   automatically.
#' @param seed Specify a random seed for the bootstrap tree.
#' @param echo Logical; should STDERR be written to the screen?
#' @param other_args Other arguments to pass to IQ tree; must be entered as a
#'   character vector with the name and value of each argument separately. For
#'   example, c("-pers", "0.2", "-nstop", "500").
#'
#' @return A single phylogenetic tree (list of class "phylo")
#'
iqtree_bs <- function(
  alignment = NULL, aln_path = NULL,
  constraint_tree, m = NULL, nt = 1, seed = 1,
  echo = FALSE, other_args = NULL) {

  # Set up temporary working directory: unique WD for each seed
  wd <- fs::path(tempdir(), glue::glue("iqtree_bs_{seed}"))
  if (fs::dir_exists(wd)) fs::dir_delete(wd)
  fs::dir_create(wd)

  if (is.null(aln_path) && is.null(alignment)) stop("Must provide either alignment or aln_path") # nolint
  if (!is.null(aln_path) && !is.null(alignment)) stop("Must provide either alignment or aln_path") # nolint

  # Write out alignment
  if (is.null(aln_path) && !is.null(alignment)) {
  assertthat::assert_that(
      is.matrix(alignment),
      msg = "alignment must be a matrix (not a list of unaligned sequences)")

  aln_path <- fs::path(wd, "alignment.phy")

  phangorn::write.phyDat(alignment, aln_path, format = "phylip")
  }

  # Write out constraint tree
  const_tree_path <- fs::path(wd, "constraint.tre")

  ape::write.tree(constraint_tree, const_tree_path)

  # Set up arguments
  iqtree_arguments <- c(
    "-s", fs::path_abs(aln_path),
    if (!is.null(nt)) "-nt", nt,
    if (!is.null(m)) "-m", m,
    if (!is.null(seed)) "-seed", seed,
    "-g", const_tree_path,
    "-bo", 1,
    "-pre", "boot",
    other_args
  )

    # Run iqtree command
  processx::run(
    "iqtree2",
    iqtree_arguments, wd = wd, echo = echo,
    # Include env variable as workaround for initial parsimony analysis
    # using all cores
    # https://github.com/iqtree/iqtree2/issues/18
    env = c("current", OMP_NUM_THREADS = "1"))

  bs_tree <- ape::read.tree(fs::path(wd, "boot.boottrees"))

  if (fs::dir_exists(wd)) fs::dir_delete(wd)

  return(bs_tree)

}

# Wrapper around read_lines that can absorb additional dummy arguments
read_lines_tar <- function(..., depends = NULL) {
  readr::read_lines(...)
}

#' Extract the best-scoring model from an IQTREE log
#'
#' @param iqtree_log Raw IQTREE log file (ending in .log) read into R
#'
#' @return Best-scoring model
#' 
extract_iqtree_mod <- function(iqtree_log) {
  iqtree_log[
    # Formatted for IQTREE2
    str_detect(iqtree_log, "Best-fit model: ")] %>%
    str_match("([^ ]+) chosen") %>%
    magrittr::extract(,2)
}

# Monophyly ----

#' Load data on Equisetum subgenera
#'
#' Filters species to those in the tree. Checks to make sure all Equisetum
#' species in the tree are included in the data.
#'
#' @param equisteum_subgen_path Path to CSV file with two columns,
#' "scientificName" and "subgenus".
#' @param sanger_tree Phylogenetic tree including Equisetum species.
#'
#' @return Tibble
#' 
load_equisetum_subgen <- function(equisteum_subgen_path, sanger_tree) {
  # Load CSV file, parse sci names
  equisetum_subgen <- read_csv(
    equisteum_subgen_path, col_types = "cc") %>%
    mutate(
      taxastand::ts_parse_names(scientificName) %>%
        select(genus = genus_name, specific_epithet),
      species = paste(genus, specific_epithet, sep = "_")) %>%
    select(species, subgenus) %>%
    unique() %>%
    # Format subgenus as it is in fossil data
    mutate(subgenus = paste3("Equisetum subgen.", subgenus))

  # Make sure all species are in tree
  tibble(
    species = sanger_tree$tip.label) %>%
    filter(str_detect(species, "Equisetum")) %>%
    left_join(equisetum_subgen, by = "species") %>%
    assert(not_na, everything()) %>%
    assert(is_uniq, species)
}

#' Add "major clade" to Sanger sampling table
#'
#' "major clade" includes pteridophyte suborder or order, and combines
#' Ophioglossales + Psilotales into one group
#'
#' @param data Tibble including columns "suborder" and "order".
#'
#' @return Tibble withh column "major_clade" added
#'
add_major_clade <- function(data) {
  data %>%
    mutate(
      major_clade = coalesce(suborder, order),
      major_clade = case_when(
        major_clade %in% c("Ophioglossales", "Psilotales") ~ "Ophioglossales + Psilotales", #nolint
        TRUE ~ major_clade
      )
    )
}

#' Make tibble summarizing sampling of Sanger dataset
#'
#' @param plastome_metadata_renamed Plastome data with final (resolved) species
#' names.
#' @param sanger_tree Sanger ML tree.
#' @param ppgi_taxonomy PPGI taxonomy
#'
#' @return Tibble with columns "species",  "genus", "order", "suborder",
#' "family"  "subfamily"  "major_clade" "outgroup"
#'
make_sanger_sampling_tbl <- function(
  plastome_metadata_renamed,
  sanger_tree, ppgi_taxonomy
) {

  # check monophyly ----
  # Make tibble of outgroup species
  og_species <-
    plastome_metadata_renamed %>%
    select(species, outgroup) %>%
    filter(outgroup == TRUE)

  # Make tibble with one row per species in Sanger sampling
  tibble(species = sanger_tree$tip.label) %>%
    # Add higher-level taxonomy
    mutate(
      genus = str_split(species, "_") %>% map_chr(1)
    ) %>%
    left_join(
      select(
        ppgi_taxonomy, order, suborder, family, subfamily, genus), by = "genus"
    ) %>%
    # Add major_clade
    add_major_clade() %>%
    # Add outgroup status
    left_join(og_species, by = "species") %>%
    mutate(outgroup = replace_na(outgroup, FALSE)) %>%
    verify(sum(outgroup) == nrow(og_species)) %>%
    # Check for match for tips with tree
    verify(all(species %in% sanger_tree$tip.label)) %>%
    verify(all(sanger_tree$tip.label %in% .$species))
}

#' Get results of monophyly test for various taxa
#'
#' @param solution Result of assessing monophyly with assess_monophy().
#' @param taxlevels Numeric vector: taxonomic levels to extract.
#'
#' @return Tibble
#'
get_result_monophy <- function(solution, taxlevels) {
  MonoPhy::GetResultMonophyly(solution, taxlevels = taxlevels) %>%
  magrittr::extract2(1) %>%
  rownames_to_column("taxon") %>%
  as_tibble() %>%
  janitor::clean_names()
}

#' Get summary of monophyly test
#'
#' @param solution Result of assessing monophyly with assess_monophy().
#' @param taxlevels Numeric vector: taxonomic levels to extract.
#'
#' @return Tibble
get_summary_monophy <- function(solution, taxlevels) {
  mp_sum <- MonoPhy::GetSummaryMonophyly(solution, taxlevels = taxlevels)

  mp_sum %>%
  magrittr::extract2(1) %>%
  rownames_to_column("var") %>%
  as_tibble() %>%
  mutate(tax_level = names(mp_sum)) %>%
  janitor::clean_names() %>%
  select(tax_level, var, taxa, tips)
}

#' Assess monophyly
#'
#' Wrapper around MonoPhy::AssessMonophyly()
#'
#' @param taxon_sampling Dataframe of taxa to assess for monophyly. Must
#' include column "species"
#' @param tree Phylogenetic tree.
#' @param og_taxa Character vector; pair of taxa to define the outgroup to
#' root the tree.
#' @param tax_levels Character vector; names of columns in `taxon_sampling`
#' to check for monophyly.
#'
#' @return List; results of MonoPhy::AssessMonophyly()
#'
assess_monophy <- function(
  taxon_sampling, tree,
  og_taxa = NULL,
  tax_levels) {
  tax_levels <- c("species", tax_levels) %>% unique()
  # Root tree
  if (!is.null(og_taxa)) {
    tree <- phytools::reroot(
      tree,
      getMRCA(tree, og_taxa)
    )
  }
  # Check monophyly
  taxon_sampling %>%
    verify("species" %in% colnames(.)) %>%
    select(species, all_of(tax_levels)) %>%
    as.data.frame() %>%
    MonoPhy::AssessMonophyly(tree, .)
}

# Dating prep ----

#' Load fossil fern data
#'
#' @param file Path to fossil fern data (CSV file)
#'
#' @return Tibble
load_fossil_data <- function(file) {
  read_csv(file) %>%
    janitor::clean_names()
}

#' Filter data on fossil calibration points
#'
#' Filters list to one point (oldest available fossil) per calibration
#' node, excludes "Incertae sedis" and known non-monophyletic taxa
#'
#' @param fossil_dates_all Tibble; fossil data read in with
#' load_fossil_data()
#'
#' @return Tibble
filter_fossil_calibration_points <- function(fossil_dates_all) {
  fossil_dates_all %>%
    # Select needed columns
    select(
      n_fos, minimum_age, node_calibrated, fossil_taxon,
      affinities_group, affinities) %>%
    # Delete any missing fossils (some records are empty that were errors)
    filter(!is.na(fossil_taxon)) %>%
    # Exclude Incertae sedis
    filter(
      str_detect(
        node_calibrated,
        regex("Incertae sedis", ignore.case = TRUE),
        negate = TRUE)
    ) %>%
    # Fix taxonomy to match pteridocat
    # Athyrium s.s. *is* Athyrium sensu pteridocat
    # Note that fossil ferns still includes Aglaomorpha, which is
    # actually a subclade of Drynaria
    mutate(
      across(c(node_calibrated, affinities),
      ~str_replace_all(., "Athyrium s.s.", "Athyrium"))
    ) %>%
    # Exclude non-monophyletic groups: Dennstaedtia, Dicksonia+Calochlaena
    filter(!affinities %in%
      c("Dennstaedtia", "Dicksonia+Calochlaena")) %>%
    # Use crown Equisetum subgen. Paramochaete (164 my) for crown Equisetum
    # since Equisetum subgen. Paramochaete is monotypic,
    # it is equivalent to dating crown Equisetum
    mutate(node_calibrated = case_when(
      node_calibrated == "crown Equisetum subgen. Paramochaete" ~
        "crown Equisetum",
      TRUE ~ node_calibrated
    )) %>%
    mutate(affinities = case_when(
      affinities == "Equisetum subgen. Paramochaete" ~
        "Equisetum",
      TRUE ~ affinities
    )) %>%
    # Keep only one oldest fossil per calibration node, no ties
    group_by(node_calibrated) %>%
    slice_max(n = 1, order_by = minimum_age, with_ties = FALSE) %>%
    ungroup()
}

#' Make a tibble mapping fossil groups (affinities) to species
#'
#' @param tree Phylogenetic tree.
#' @param fossil_calibration_points Tibble of fossil calibration points.
#' @param ppgi_taxonomy Pteridophyte phylogeny group I taxonomy.
#' @param equisetum_subgen Subgenera of Equisetum and their species.
#' @param plastome_metadata_renamed Metada for plastome sequences, including
#' species and outgroup status.
#' @param include_algaomorpha Logical; should the Aglaomorpha subclade of
#' Drynaria be included?
#'
#' @return Tibble with two columns, "affinities" and "species"
#'
make_fossil_species_map <- function(
  tree, fossil_calibration_points, ppgi_taxonomy,
  equisetum_subgen, plastome_metadata_renamed,
  include_algaomorpha = TRUE) {

  # Make tibble of tips with genus and species
  tip_tbl <- tibble(species = tree$tip.label) %>%
    mutate(genus = str_split(species, "_") %>% map_chr(1)) %>%
    # Modify for Polypodium s.l.: includes Pleurosoriopsis
    mutate(genus = case_when(
      genus == "Polypodium" ~ "Polypodium s.l.",
      genus == "Pleurosoriopsis" ~ "Polypodium s.l.",
      TRUE ~ genus
    )) %>%
    assert(is_uniq, species) %>%
    assert(not_na, everything())

  # Filter PPGI taxonomy to only genera in the tree
  ppgi_taxonomy_in_tree <-
    ppgi_taxonomy %>%
    # Modify for Polypodium s.l.: includes Pleurosoriopsis
    mutate(genus = case_when(
      genus == "Polypodium" ~ "Polypodium s.l.",
      genus == "Pleurosoriopsis" ~ "Polypodium s.l.",
      TRUE ~ genus
    )) %>%
    inner_join(unique(select(tip_tbl, genus)), by = "genus") %>%
    unique() %>%
    assert(is_uniq, genus)

  # Filter equisetum subgenera to only species in the tree
  equisetum_subgen_in_tree <-
    equisetum_subgen %>%
    inner_join(select(tip_tbl, species), by = "species")

  # Make tibble of deeper groups (Euphyllophytes and Tracheophytes)
  # Needs to include outgroup taxa
  bryo_genera <- c("Anthoceros", "Physcomitrium", "Marchantia")
  lyco_genera <- c("Isoetes", "Lycopodium", "Selaginella")

  og_deep_clades <-
    plastome_metadata_renamed %>%
    filter(outgroup == TRUE) %>%
    select(species) %>%
    mutate(
      clade_1 = if_else(
        str_detect(species, paste(bryo_genera, collapse = "|"), negate = TRUE),
        "Tracheophytes", NA_character_),
      clade_2 = if_else(
        str_detect(species, 
          paste(c(bryo_genera, lyco_genera), collapse = "|"), negate = TRUE),
        "Euphyllophytes", NA_character_)
    ) %>%
    select(species, contains("clade")) %>%
    arrange(clade_2, clade_1, species)

  deep_clades <-
  tip_tbl %>%
    anti_join(og_deep_clades, by = "species") %>%
    select(-genus) %>%
    mutate(clade_1 = "Tracheophytes", clade_2 = "Euphyllophytes") %>%
    bind_rows(og_deep_clades) %>%
    arrange(clade_2, clade_1, species) %>%
    verify(all(species %in% tip_tbl$species)) %>%
    verify(all(tip_tbl$species %in% species)) %>%
    verify(nrow(.) == nrow(tip_tbl))

  # Make tibble of Aglaomorpha (subclade of Drynaria) to include
  if (include_algaomorpha == TRUE) {
    aglaomorpha_tbl <-
    ape::getMRCA(
      tree, c("Drynaria_meyeniana", "Drynaria_novoguineensis")) %>%
      get_children(tree, .) %>%
      tibble(
        species = .,
        genus = "Aglaomorpha"
      ) %>% 
    # Make sure Aglaomorpha is formatted as expected:
    # check that immediate sister species
    # *outside* Aglaomorpha are *not* included
    verify(!"Drynaria_mollis" %in% .$species) %>%
    verify(!"Drynaria_fortunei" %in% .$species)
  
    tip_tbl <- bind_rows(tip_tbl, aglaomorpha_tbl)
  }
    
  fossil_calibration_points %>%
    select(affinities) %>%
    # Split affinities that are composed of multiple taxa separated by '+'
    mutate(aff_split = affinities) %>%
    separate_rows(aff_split, sep = "\\+") %>%
    # Affinities comprise subgenus (Equisetum only), subfamily, family, order
    # - join genus by subfamily
    left_join(
      select(
        ppgi_taxonomy_in_tree,
        aff_split = subfamily, genus_1 = genus),
      by = "aff_split"
    ) %>%
    # - join genus by family
    left_join(
      select(
        ppgi_taxonomy_in_tree,
        aff_split = family, genus_2 = genus),
      by = "aff_split"
    ) %>%
    # - join genus by order
    left_join(
      select(
        ppgi_taxonomy_in_tree,
        aff_split = order, genus_3 = genus),
      by = "aff_split"
    ) %>%
    mutate(genus = coalesce(genus_1, genus_2, genus_3, aff_split)) %>%
    select(affinities, aff_split, genus) %>%
    unique() %>%
    assert(not_na, everything()) %>%
    # - join species by genus
    left_join(tip_tbl, by = "genus") %>%
    # Join equisetum species by subgenus
    left_join(
      select(
        equisetum_subgen_in_tree,
        species_2 = species,
        aff_split = subgenus
      ),
      by = "aff_split"
    ) %>%
    # Join deeper groups by species
    left_join(
      select(deep_clades, species_3 = species, aff_split = clade_1),
      by = "aff_split"
    ) %>%
    left_join(
      select(deep_clades, species_4 = species, aff_split = clade_2),
      by = "aff_split"
    ) %>%
    mutate(species = coalesce(species, species_2, species_3, species_4)) %>%
    select(affinities, species) %>%
    unique() %>%
    assert(not_na, everything())
}

#' Parse fossil calibration data in Sundue and Testo 2016 SI
#'
#' Testo WL, Sundue MA (2016) A 4000-species dataset provides new insight into
#' the evolution of ferns. Molecular Phylogenetics and Evolution 105:200–211.
#' https://doi.org/10.1016/j.ympev.2016.09.003
#'
#' @param testo_sundue_2016_si_path Path to Sundue and Testo 2016 SI file
#' (xlsx format).
#'
#' @return Tibble with fossil calibration dates
#' 
parse_ts_calibrations <- function(testo_sundue_2016_si_path) {
readxl::read_excel(
  testo_sundue_2016_si_path,
  # Divergence dates are in the 5th sheet
  sheet = 3
) %>%
  janitor::clean_names() %>%
  mutate(
    stem_crown = str_to_lower(stem_crown),
    clade = 
      str_replace_all(clade, "Aglaomorpha heraclea", "Drynaria heraclea") %>%
      str_replace_all("Alsophila/Cyathea clade", "Alsophila+Cyathea") %>%
      str_replace_all("Isoetales", "Isoëtales") %>%
      str_replace_all(" ", "_")) %>%
  transmute(
    minimum_age = age,
    node_calibrated = paste(stem_crown, clade),
    fossil_taxon = fossil,
    affinities_group = stem_crown,
    affinities = clade
  ) %>%
  mutate(across(
    c(node_calibrated, affinities),
    ~str_replace_all(., "Polypodium", "Polypodium s.l.")))
}

#' Make a tibble mapping fossil groups (affinities) to species
#' for Testo and Sundue 2016 data
#'
#' @param tree Phylogenetic tree.
#' @param fossil_calibration_points Tibble of fossil calibration points.
#' @param ppgi_taxonomy Pteridophyte phylogeny group I taxonomy.
#' @param plastome_metadata_renamed Metada for plastome sequences, including
#' species and outgroup status.
#'
#' @return Tibble with two columns, "affinities" and "species"
#'
make_ts_fossil_species_map <- function(
  tree, fossil_calibration_points, ppgi_taxonomy,
  plastome_metadata_renamed) {

  # Make tibble of tips with genus and species
  tip_tbl <- tibble(species = tree$tip.label) %>%
    mutate(genus = str_split(species, "_") %>% map_chr(1)) %>%
    # Modify for Polypodium s.l.: includes Pleurosoriopsis
    mutate(genus = case_when(
      genus == "Polypodium" ~ "Polypodium s.l.",
      genus == "Pleurosoriopsis" ~ "Polypodium s.l.",
      TRUE ~ genus
    )) %>%
    assert(is_uniq, species) %>%
    assert(not_na, everything())

  # Filter PPGI taxonomy to only genera in the tree
  ppgi_taxonomy_in_tree <-
    ppgi_taxonomy %>%
    # Modify for Polypodium s.l.: includes Pleurosoriopsis
    mutate(genus = case_when(
      genus == "Polypodium" ~ "Polypodium s.l.",
      genus == "Pleurosoriopsis" ~ "Polypodium s.l.",
      TRUE ~ genus
    )) %>%
    inner_join(unique(select(tip_tbl, genus)), by = "genus") %>%
    unique() %>%
    assert(is_uniq, genus)

  # Make tibble of deeper clade (Euphyllophytes)
  # Needs to include outgroup taxa
  bryo_genera <- c("Anthoceros", "Physcomitrium", "Marchantia")
  lyco_genera <- c("Isoetes", "Lycopodium", "Selaginella")

  og_deep_clades <-
    plastome_metadata_renamed %>%
    filter(outgroup == TRUE) %>%
    select(species) %>%
    mutate(
      clade = if_else(
        str_detect(species, 
          paste(c(bryo_genera, lyco_genera), collapse = "|"), negate = TRUE),
        "Euphyllophytes", NA_character_)
    ) %>%
    select(species, contains("clade")) %>%
    arrange(clade, species)

  deep_clades <-
  tip_tbl %>%
    anti_join(og_deep_clades, by = "species") %>%
    select(-genus) %>%
    mutate(clade = "Euphyllophytes") %>%
    bind_rows(og_deep_clades) %>%
    verify(all(species %in% tip_tbl$species)) %>%
    verify(all(tip_tbl$species %in% species)) %>%
    verify(nrow(.) == nrow(tip_tbl))

  fossil_calibration_points %>%
    select(affinities) %>%
    # Split affinities that are composed of multiple taxa separated by '+'
    mutate(aff_split = affinities) %>%
    separate_rows(aff_split, sep = "\\+") %>%
    # Affinities comprse family, order, some species
    # - join genus by family
    left_join(
      select(
        ppgi_taxonomy_in_tree,
        aff_split = family, genus_1 = genus),
      by = "aff_split"
    ) %>%
    # - join genus by order
    left_join(
      select(
        ppgi_taxonomy_in_tree,
        aff_split = order, genus_2 = genus),
      by = "aff_split"
    ) %>%
    # - join genus by genus
    left_join(
      select(
        ppgi_taxonomy_in_tree,
        aff_split = genus, genus_3 = genus),
      by = "aff_split"
    ) %>%
    mutate(genus = coalesce(genus_3, genus_2, genus_1)) %>%
    select(affinities, aff_split, genus) %>%
    unique() %>%
    # - join species by genus
    left_join(tip_tbl, by = "genus") %>%
    # Join deeper groups by species
    left_join(
      select(deep_clades, species_2 = species, aff_split = clade),
      by = "aff_split"
    ) %>%
    # Fill in affinities at species level
    left_join(
      select(tip_tbl, aff_split = species, species_3 = species),
      by = "aff_split"
    ) %>%
    mutate(species = coalesce(species, species_2, species_3)) %>%
    select(affinities, species) %>%
    unique() %>%
    assert(not_na, everything())
}

#' Test monophyly of fossil groups
#'
#' @param tree Phylogenetic tree, should be rooted.
#' @param fossil_node_species_map Tibble with two columns,
#' "affinities" and "species".
#' @param fossil_affinity_select Name of fossil group to test.
#'
#' @return Tibble with monophyletic status of selected fossil group
#' 
mono_test_fossil <- function(
  tree,
  fossil_node_species_map,
  fossil_affinity_select
) {
  # Make tibble with taxon sampling of fossil groups
  taxon_sampling <-
    tibble(species = tree$tip.label) %>%
    left_join(
      filter(fossil_node_species_map, affinities == fossil_affinity_select),
      by = "species"
    )
  # Test monophyly of selected group
  assess_monophy(
    taxon_sampling = taxon_sampling,
    tree = tree,
    tax_levels = "affinities"
  ) %>%
  get_result_monophy(., 1)
}

#' Get a pair of tips that define a clade in a phylogenetic tree
#'
#' @param tree Phylogenetic tree, must be rooted.
#' @param node Number of a node in the tree.
#'
#' @return Character vector; a pair of tips whose MRCA is `node`
#' 
get_spanning_tips <- function(tree, node) {
  # Tree must be rooted
  assertthat::assert_that(ape::is.rooted(tree))
  # Ladderize tree
  tree <- ape::ladderize(tree)
  # Get vector of tips in ladderized order
  all_tips_ladder_ord <- get_tips_in_ape_plot_order(tree)
  # Get spanning tips, not in ladderized order
  spanning_tips_unord <- phangorn::Descendants(tree, node, "tips") %>%
    magrittr::extract2(1) %>%
    magrittr::extract(tree$tip.label, .)
  # Put 'final' spanning tips in ladderized order
  final_tips <-
    all_tips_ladder_ord[all_tips_ladder_ord %in% spanning_tips_unord]
  # Return first and last spanning tips in ladderized order
  if (length(final_tips) > 2) return(
    c(dplyr::first(final_tips), dplyr::last(final_tips)))
  final_tips
}

#' Count the number of tips descending from a node defined by a pair of tips 
#' in a phylogenetic tree
#'
#' @param tree Phylogenetic tree.
#' @param tips Character vector; a pair of tips in the tree.
#'
#' @return Number of terminal tips in the clade that is defined by the MRCA of 
#' `tips`
#' 
check_num_tips <- function(tree, tips) {
  if (is.null(tips)) return(NA)
  getMRCA(tree, tips) %>%
    phangorn::Descendants(tree, node = ., "tips") %>%
    magrittr::extract2(1) %>%
    length()
}

#' Get parent node from a tip taxon
#'
#' @param tree Phylogenetic tree.
#' @param species Single tip of the tree.
#'
#' @return Number of the parent node of the species
#' 
get_parent <- function(tree, species) {
  node <- which(tree$tip.label == species)
  phangorn::Ancestors(tree, node, type = "parent")
}

#' Get child tips (taxa) from parent node
#'
#' @param tree Phylogenetic tree.
#' @param node Number of a node in the tree.
#'
#' @return Character vector; tips that descend from `node`
#' 
get_children <- function(tree, node) {
  # Tree must be rooted
  assertthat::assert_that(ape::is.rooted(tree))
  # Get spanning tips, not in ladderized order
  phangorn::Descendants(tree, node, "tips") %>%
    magrittr::extract2(1) %>%
    magrittr::extract(tree$tip.label, .)
}

#' Make tibble of manual spanning tips
#' 
#' This should always be done after inspecting the actual tree used for dating
#' 
#' @param data_set Choose data set; either "this_study" or "ts2016"
#' (Testo and Sundue 2016 fossil calibration points)
#' @return Tibble with columns "affinities", "tip_1_manual", and "tip_2_manual".
#' - "affinities" is the name of the group that corresponds to
#' a fossil calibration point.
#' - The MRCA of "tip_1_manual" and "tip_2_manual" define each group
#' in "affinities"
define_manual_spanning_tips <- function(data_set = c("this_study", "ts2016")) {
  switch(
    data_set,
    "this_study" = tribble(
      ~affinities, ~affinities_group, ~tip_1_manual, ~tip_2_manual,
      "Pleopeltis", "crown", "Pleopeltis_bombycina", "Pleopeltis_conzattii",
      "Polypodium s.l.", "stem", "Pleurosoriopsis_makinoi", "Polypodium_pellucidum",
      "Cyathea+Alsophila+Gymnosphaera", "crown", "Alsophila_poolii", "Cyathea_epaleata"
      ),
    "ts2016" = tribble(
      ~affinities, ~affinities_group, ~tip_1_manual, ~tip_2_manual,
      "Alsophila+Cyathea", "stem", "Cyathea_minuta", "Alsophila_capensis",
      "Diplazium+Athyrium", "stem", "Ephemeropteris_tejeroi", "Diplazium_caudatum",
      "Pleopeltis", "crown", "Pleopeltis_bombycina", "Pleopeltis_conzattii",
      "Polypodium s.l.", "stem", "Pleurosoriopsis_makinoi", "Polypodium_pellucidum"),
    stop("Must choose either 'this_study' or 'ts2016'")
  )
}

#' Get tips that define clades corresponding to each fossil
#' calibration point
#' 
#' Filters out redundant calibrations (fossils calibrating the same node)
#'
#' @param fossil_node_species_map Tibble in long format with two columns,
#' "affinities" and "species".
#' @param sanger_tree_rooted Rooted phylogeny.
#' @param fossil_calibration_points Fossil calibration points read in with
#' load_fossil_calibration_points().
#' @param manual_spanning_tips Tibble of tips that define clades corresponding
#' to fossils prepared by hand; used for groups that are non-monophyletic in
#' sanger_tree_rooted.
#'
#' @return Tibble of fossil calibration points with two columns added, "tip_1"
#' and "tip_2" that define the clades corresponding to each fossil
#' calibration point
#'
get_fossil_calibration_tips <- function(
  fossil_node_species_map,
  sanger_tree_rooted,
  fossil_calibration_points,
  manual_spanning_tips
) {
  # Make sure phylogeny is rooted
  assertthat::assert_that(ape::is.rooted(sanger_tree_rooted))

  # Check monophyly of each group
  fossil_node_monophy <-
  map_df(
    sort(unique(fossil_node_species_map$affinities)),
    ~mono_test_fossil(
      tree = sanger_tree_rooted,
      fossil_node_species_map = fossil_node_species_map,
      fossil_affinity_select = .
    )
  ) %>%
    mutate(across(matches("mrca|number|delta"), parse_number)) %>%
    rename(affinities = taxon) %>%
    # Set MRCA to NA if not monophyletic
    # (so a manual MRCA can be added later)
    mutate(mrca = case_when(
      monophyly == "No" ~ NaN,
      TRUE ~ mrca
    ))

  # Double check that all non-monophyletic taxa are specified in
  # manual_spanning_tips (both affinities and affinities_group must match)
  aff_non_mono <-
  fossil_node_monophy %>%
    filter(monophyly == "No") %>%
    select(affinities) %>%
    unique()

  fossil_calibration_points %>%
    inner_join(aff_non_mono, by = "affinities") %>%
    select(affinities, affinities_group) %>%
    anti_join(
      manual_spanning_tips, by = c("affinities", "affinities_group")) %>%
    verify(
      nrow(.) == 0,
      success_fun = success_logical,
      error_fun = err_msg("Not exact match between non-monophyletic groups and manual spanning tips") # nolint
    )

  # Make tibble of stem MRCA for monotypic calibration groups
  monotypic_stem_mrca_tib <-
  fossil_node_monophy %>%
    filter(monophyly == "Monotypic") %>%
    mutate(affinities_group = "stem") %>%
    left_join(fossil_node_species_map, by = "affinities") %>%
    select(affinities, affinities_group, species) %>%
    mutate(
      monotypic_stem_mrca = map_dbl(species, ~get_parent(sanger_tree_rooted, .))
    ) %>%
    # Check that species is amongst descendents from MRCA
    mutate(
      children = map(monotypic_stem_mrca, ~get_children(sanger_tree_rooted, .)),
      sp_in_children = map2_lgl(species, children, ~magrittr::is_in(.x, .y))
    ) %>%
    assert(isTRUE, sp_in_children) %>%
    select(affinities, affinities_group, monotypic_stem_mrca)

  # Calculate mrca for manually specified groups (non-monophyletic taxa)
  manual_mrca <-
    manual_spanning_tips %>%
    mutate(
      mrca_manual = map2_dbl(
        tip_1_manual, tip_2_manual,
        ~ape::getMRCA(sanger_tree_rooted, c(.x, .y))
        ),
      stem_mrca_manual = phangorn::Ancestors(
        sanger_tree_rooted, mrca_manual, type = "parent")
    ) %>%
    select(-tip_1_manual, -tip_2_manual)

  # Make tibble of tips spanning each fossil group
  # in "long" format with spanning tips as list-col
  spanning_tips_long <-
    # Start with fossil calibration points
    fossil_calibration_points %>%
    # Add MRCA and monophyly status
    left_join(
      select(fossil_node_monophy, affinities, monophyly, number_tips, mrca),
      by = "affinities"
    ) %>%
    # Add 'stem MRCA': the parent of each MRCA, used for stem groups
    mutate(
      stem_mrca = map_dbl(
        mrca, ~phangorn::Ancestors(sanger_tree_rooted, ., "parent"))
    ) %>%
    # Add stem MRCA for monotypic stem groups
    # (these will still lack "normal" mrca)
    left_join(
      monotypic_stem_mrca_tib,
      by = c("affinities", "affinities_group")
    ) %>%
    mutate(
      stem_mrca = coalesce(stem_mrca, monotypic_stem_mrca)
    ) %>%
    select(-monotypic_stem_mrca) %>%
    # Add MRCA and stem MRCA for manually specified groups
    left_join(
      manual_mrca,
      by = c("affinities", "affinities_group")
    ) %>%
    mutate(
      mrca = coalesce(mrca, mrca_manual),
      stem_mrca = coalesce(stem_mrca, stem_mrca_manual)
    ) %>%
    select(-mrca_manual, -stem_mrca_manual) %>%
    # monotypic groups have NA for mrca, but not stem_mrca
    assert(not_na, stem_mrca) %>%
    mutate(
      # Add tips spanning each crown group for non-monotypic taxa
      rep_tips_crown = case_when(
        monophyly != "Monotypic" ~ map(
          mrca, ~get_spanning_tips(sanger_tree_rooted, .))
      ),
      # Add tips spanning each stem group for all taxa
      rep_tips_stem = map(stem_mrca, ~get_spanning_tips(sanger_tree_rooted, .))
    ) %>%
    # Add double check on number of tips descendend from crown group
    # spanning tips
    mutate(
      num_tips_check = map_dbl(
        rep_tips_crown, ~check_num_tips(sanger_tree_rooted, .))
    ) %>%
    # Use stem or crown tips as appropriate
    mutate(rep_tips = case_when(
      affinities_group == "crown" ~ rep_tips_crown,
      affinities_group == "stem" ~ rep_tips_stem,
    )) %>%
    select(-rep_tips_crown, -rep_tips_stem)

  # Check that spanning tips cover all expected species
  # or no MRCA exists in case of monotypic groups
  spanning_tips_long %>%
    filter(monophyly %in% c("Yes", "Monotypic")) %>%
    verify(
      all(number_tips == num_tips_check | is.na(mrca)),
      error_fun = err_msg("Spanning tips do not match fossil group"),
      success_fun = success_logical)

  # Convert rep tips from list-col to two columns, "tip_1" and "tip_2"
  # Also drop any remaining redundant calibration points
  spanning_tips_nulls_gone <-
    spanning_tips_long %>%
    select(-num_tips_check) %>%
    # drops NULLs so will need to rejoin later
    unnest(rep_tips) %>%
    group_by(node_calibrated) %>%
    mutate(n_tip = 1:n() %>%
      paste0("tip_", .)) %>%
    ungroup() %>%
    pivot_wider(values_from = rep_tips, names_from = n_tip) %>%
    assert(is_uniq, node_calibrated) %>%
    assert(not_na, node_calibrated) %>%
    # Filter out redundant calibration points (same tips)
    assert(not_na, tip_1, tip_2) %>%
    group_by(tip_1, tip_2) %>%
    slice_max(n = 1, order_by = minimum_age, with_ties = FALSE) %>%
    ungroup()

  # Add back in monotypic groups (rep_tips are NULL)
  # to get to final "spanning_tips" tbl
  spanning_tips <-
    # Combine spanning tips with rep_tips NULL to widened data
    spanning_tips_long %>%
    filter(map_lgl(rep_tips, is.null)) %>%
    select(-rep_tips) %>%
    bind_rows(spanning_tips_nulls_gone) %>%
    assert(is_uniq, node_calibrated) %>%
    assert(not_na, node_calibrated) %>%
    select(-num_tips_check)

  ## Check results ##
  # Double check that there are no redundant tip sets
  spanning_tips %>%
    filter(!is.na(tip_1)) %>%
    add_count(tip_1, tip_2) %>%
    verify(
      all(n == 1), success_fun = success_logical,
      error_fun = err_msg("Redundant spanning tips detected"))

  # Make sure manual tips cover all non-monophyletic groups
  spanning_tips %>%
    filter(monophyly == "No") %>%
    anti_join(manual_spanning_tips, by = "affinities") %>%
    verify(nrow(.) == 0, success_fun = success_logical,
    error_fun = err_msg("Manually specified tips do not cover all non-monophyletic groups")) #nolint

  # Check that stem age is older than crown age for each affinity with
  # both crown and stem
  aff_with_stem_and_crown <-
  spanning_tips %>%
    group_by(affinities) %>%
    add_count() %>%
    ungroup() %>%
    filter(n > 1)
  
  if (nrow(aff_with_stem_and_crown) > 0) {
    aff_with_stem_and_crown %>%
    select(minimum_age, affinities_group, affinities) %>%
    pivot_wider(names_from = affinities_group, values_from = minimum_age) %>%
    mutate(stem_older = stem > crown) %>%
    verify(
      all(stem_older == TRUE),
      success_fun = success_logical,
      error_fun = err_msg("At least one crown age is older than stem age"))
  }

  # Run final checks and output results
  spanning_tips %>%
    assert(not_na,
      minimum_age, node_calibrated, fossil_taxon, affinities_group,
      affinities, monophyly, number_tips, tip_1, tip_2) %>%
    assert(is_uniq, node_calibrated)
}

#' Make tibble of times for calibrating the root of a tree when
#' dating with treePL
#'
#' @param tree Phylogenetic tree, must be rooted.
#' @param node_name Name to use for the root node.
#' @param time Time to use for calibrating the root node.
#' @param tip_1 Tip of the tree used to define clade including all tips in tree.
#' @param tip_2 Tip of the tree used to define clade including all tips in tree.
#'
#' @return Tibble with columns "mrca", "min", and "max" formatted like
#' lines "mrca", "min", and "max" in a treePL config file
#'
calibrate_root_node <- function(tree, node_name, time, tip_1, tip_2) {
  # Tree must be rooted
  assertthat::assert_that(ape::is.rooted(tree))
  # Check that the clade defined by tip_1 and tip_2 includes all tips
  # in the tree
  tips_descended_from_root <-
    ape::getMRCA(tree, c(tip_1, tip_2)) %>%
      get_children(tree, .)

  assertthat::assert_that(
    length(tips_descended_from_root) == ape::Ntip(tree),
    msg = "Tips used to define root node do not include MRCA of all tips in tree") #nolint

  tribble(
    ~mrca, ~min, ~max, ~taxon_1, ~taxon_2,
    glue::glue("mrca = {node_name} {tip_1} {tip_2}"),
    glue::glue("min = {node_name} {time}"),
    glue::glue("max = {node_name} {time}"),
    tip_1,
    tip_2
  ) %>%
  mutate(across(everything(), as.character))
}

#' Format fossil calibration points so they can be used for treePL
#' 
#' See https://github.com/blackrim/treePL/wiki/Quick-run
#'
#' @param fossil_calibration_tips Tibble of fossil calibration points.
#' @param root_tibble Tibble with calibration data for root node.
#'
#' @return Tibble with columns "mrca", "min", and "max" formatted like
#' lines "mrca", "min", and "max" in a treePL config file
#' 
format_calibrations_for_treepl <- function(
  fossil_calibration_tips, root_tibble) {

  # Helper to check for presence of non-alphabetic or underscore character in a
  # string
  all_alpha_or_underscore <- function(x) {
    stringr::str_detect(x, "[^_a-zA-Z]", negate = TRUE)
    }

  fossil_calibration_tips %>%
    # Modify name of calibrated nodes to only use alphabet or underscore
    mutate(
      node_calibrated = str_replace_all(node_calibrated, " ", "_") %>%
        str_replace_all("\\+", "_") %>%
        str_remove_all("\\.") %>%
        # Tranform non-ascii characters
        stringi::stri_trans_general("latin-ascii")
    ) %>%
    assert(is_uniq, node_calibrated) %>%
    assert(all_alpha_or_underscore, node_calibrated) %>%
    transmute(
      mrca = glue::glue("mrca = {node_calibrated} {tip_1} {tip_2}") %>%
        as.character(),
      min = glue::glue("min = {node_calibrated} {minimum_age}") %>%
        as.character(),
      taxon_1 = tip_1,
      taxon_2 = tip_2
    ) %>%
    bind_rows(root_tibble) %>%
    assert(not_na, mrca, taxon_1, taxon_2) %>%
    assert(is_uniq, mrca) %>%
    select(mrca, min, max, taxon_1, taxon_2)
}

# Etc ----
# This function can be called inside of other functions to check
# if the names of the input match the names of the arguments
# Need to provide match.call() as the input
# check_args(match.call())
check_args <- function(call_match) {
  call_names <- as.character(call_match)
  arg_names <- names(as.list(call_match))
  stopifnot(
    "Names of input must match names of arguments" = isTRUE(all.equal(call_names[-1], arg_names[-1]))
  )
}

# Generate custom error message for assertr
# https://github.com/ropensci/assertr/issues/70
err_msg <- function(msg) stop(msg, call. = FALSE)

# Write fasta and return path
write_fasta_tar <- function(x, file, ...) {
  ape::write.FASTA(x = x, file = file, ...)
  file
  }

# Write phylogeny and return path
write_tree_tar <- function(phy, file, ...) {
  ape::write.tree(phy = phy, file = file, ...)
  file
  }

# Write CSV and return path
write_csv_tar <- function(x, file, ...) {
  readr::write_csv(x = x, file = file, ...)
  file
}

#' Calculate base frequences *per sample* in an alignment
#'
#' @param seqs DNA list (list of class "DNAbin")
#'
#' @return Number of ACTG bases in seqs
count_actg <- function(seqs) {
  count_all <- ape::base.freq(seqs, freq = TRUE, all = TRUE)
  sum(count_all[names(count_all) %in% c("a", "c", "t", "g")])
}

#' Download sequences from GenBank by accession number
#'
#' Wrapper for ape::read.GenBank() that can handle large numbers of accessions
#'
#' @param accessions Vector of GenBank accessions (e.g., "AY178864")
#' @param return_seqtbl Logical; should the results be returned as a DNA-sequence tibble (seqtbl)?
#' If FALSE, will return as list of class 'DNAbin'.
#'
#' @return Tibble (seqtbl) or list of class DNAbin
read_genbank <- function(accessions, return_seqtbl = TRUE) {

   # Split accessions into chunks
  chunk_size <- 50
  n <- length(accessions)
  r <- rep(1:ceiling(n/chunk_size), each=chunk_size)[1:n]
  accessions_list <- split(accessions, r) %>% 
    magrittr::set_names(NULL)

  # Download sequences in each chunk
  seqs_list <- purrr::map(accessions_list, ~ape::read.GenBank(., species.names = FALSE))

  # Combine results
  results <- do.call(c, seqs_list)

  # Optionally convert to seqtbl
  if(isTRUE(return_seqtbl)) results <- dnabin_to_seqtbl(results)

  results

}

# Vectorized check for NULL in a list
not_null <- function(x) {
  assertthat::assert_that(is.list(x))
  !purrr::map_lgl(x, is.null)
}

# Convert ape::DNAbin format to Biostrings::DNAStringSet format,
# optionally removing gaps
DNAbin_to_DNAstringset <- function (seqs, remove_gaps = TRUE) {
  if(isTRUE(remove_gaps)) {
  seqs %>% as.list() %>% as.character %>% 
      lapply(.,paste0,collapse="") %>% 
      lapply( function (x) gsub("-", "", x)) %>% 
      unlist %>% Biostrings::DNAStringSet()
  } else {
    seqs %>% as.list() %>% as.character %>% 
      lapply(.,paste0,collapse="") %>% 
      unlist %>% Biostrings::DNAStringSet()
  }
}

#' paste3
#'
#' Paste while removing NAs
#'
#' Removes NAs from pasted elements, but if ALL elements are NA, the result is NA.
#'
#' Shamelessly copied from
#' \url{https://stackoverflow.com/questions/13673894/suppress-nas-in-paste}
#' @param ... Strings to paste
#' @param sep Character used to separate pasted strings
#' @examples
#' paste3(c("a", "b", "c", NA), c("A","B", NA, NA))
#' @export
paste3 <- function(..., sep=" ") {
  L <- list(...)
  L <- lapply(L,function(x) {x[is.na(x)] <- ""; x})
  ret <-gsub(paste0("(^",sep,"|",sep,"$)"),"",
             gsub(paste0(sep,sep),sep,
                  do.call(paste,c(L,list(sep=sep)))))
  is.na(ret) <- ret==""
  ret
}

# Quiet version of rgnparser::gn_parse_tidy()
gn_parse_tidy_quiet <- function(...) {
  suppressMessages(
    rgnparser::gn_parse_tidy(...)
  )
}

#' Get tips of a phylogenetic tree in their plotted order
#'
#' After re-rooting a tree, the order of tips when the tree
#' is plotted no longer match the order of $tip.label. Use
#' this function to get tips in the order they are plotted.
#' @param tree List of class "phylo"
#' @return Character vector
get_tips_in_ape_plot_order <- function (tree) {
  assertthat::assert_that(inherits(tree, "phylo"))
  # First filter out internal nodes
  # from the the second column of the edge matrix
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  # Use this vector to extract the tips in the right order
  tree$tip.label[ordered_tips]
}
