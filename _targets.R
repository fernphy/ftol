library(targets)
library(tarchetypes)
source("R/packages.R")
source("R/functions.R")

# Specify path to folder with raw data
data_raw <- "/data_raw"

tar_plan(
  
  # Load data ----
  # Data for resolving taxonomic names:
  # - Catalog of Life database
  tar_file(col_data_path, path(data_raw, "2021-08-25_dwca/Taxon.tsv")),
  col_data = load_col(col_data_path),
  # - World Ferns taxonomic data
  world_ferns_data = extract_fow_from_col(col_data),
  # - Modified PPGI taxonomy
  # with new genera and slightly different treatments following World Ferns list
  tar_file(ppgi_taxonomy_path, path(data_raw, "ppgi_taxonomy_mod.csv")),
  ppgi_taxonomy = read_csv(ppgi_taxonomy_path),
  # List of plastid coding genes from Wei et al 2017
  tar_file(target_sanger_genes_path, path(data_raw, "wei_2017_coding_genes.txt")),
  target_sanger_genes = read_lines(target_sanger_genes_path),
  # Outgroup plastome accessions
  tar_file(plastome_outgroups_path, path(data_raw, "plastome_outgroups.csv")),
  plastome_outgroups = read_csv(plastome_outgroups_path),
  # Calibration dates after Testo and Sundue 2016
  tar_file(plastome_calibration_dates_path, path(data_raw, "testo_sundue_2016_calibrations.csv")),
  plastid_calibration_dates = load_calibration_dates(plastome_calibration_dates_path),
  # Manually selected synonyms for resolving names of plastid genes
  tar_file(sanger_names_with_mult_syns_select_path, path(data_raw, "genbank_names_with_mult_syns_select.csv")),
  sanger_names_with_mult_syns_select = read_csv(sanger_names_with_mult_syns_select_path),
  
  # Download individual plastid sequences from GenBank (Sanger sequences) ----
  # Define variables used in plan
  # - Target plastid fern genes to download
  target_genes = c("atpA", "atpB", "rbcL", "rps4"),
  # - Most recent date cutoff for sampling genes
  date_cutoff = "2021/09/01",
  # Download fern plastid gene fasta files
  tar_target(
    raw_fasta,
    fetch_fern_gene(target_genes, end_date = date_cutoff),
    pattern = map(target_genes),
    deployment = "main" # don't run in parallel, or will get HTTP status 429 errors
  ),
  # Download fern plastid gene metadata
  tar_target(
    raw_meta,
    fetch_fern_metadata(target_genes, end_date = date_cutoff),
    pattern = map(target_genes),
    deployment = "main"
  ),
  
  # Taxonomic name resolution ----
  # Download species names from NCBI
  ncbi_names_raw = raw_meta %>% pull(taxid) %>% unique %>% fetch_taxonomy,
  # Clean NCBI species names
  ncbi_names_full = clean_ncbi_names(ncbi_names_raw),
  # Exclude invalid names (hybrids, taxa not identified to species level)
  ncbi_names_query = exclude_invalid_ncbi_names(ncbi_names_full),
  # Parse reference names
  wf_ref_names = tt_parse_names(unique(world_ferns_data$scientific_name)),
  # Resolve names, round 1: NCBI accepted scientific names
  ncbi_names_query_round_1 = select_ncbi_names_round_1(ncbi_names_query),
  # - match names to reference
  match_results_raw_round_1 = tt_match_names(
    query = ncbi_names_query_round_1$scientific_name, 
    reference = wf_ref_names,
    max_dist = 5, match_no_auth = TRUE, match_canon = TRUE) %>%
    as_tibble(),
  # - classify matching results
  match_results_classified_round_1 = tt_classify_result(match_results_raw_round_1),
  # - resolve synonyms
  match_results_resolved_round_1 = tt_resolve_synonyms(match_results_classified_round_1, world_ferns_data),
  # Resolve names, round 2: NCBI synonym scientific names
  ncbi_names_query_round_2 = select_ncbi_names_round_2(match_results_resolved_round_1, ncbi_names_query),
  match_results_raw_round_2 = tt_match_names(
    query = ncbi_names_query_round_2$scientific_name, 
    reference = wf_ref_names,
    max_dist = 5, match_no_auth = TRUE, match_canon = TRUE) %>%
    as_tibble(),
  match_results_classified_round_2 = tt_classify_result(match_results_raw_round_2),
  match_results_resolved_round_2 = tt_resolve_synonyms(match_results_classified_round_2, world_ferns_data),
  # Resolve names, round 3: NCBI species without author
  ncbi_names_query_round_3 = select_ncbi_names_round_3(match_results_resolved_round_1, match_results_resolved_round_2, ncbi_names_query),
  match_results_raw_round_3 = tt_match_names(
    query = ncbi_names_query_round_3$species, 
    reference = wf_ref_names,
    max_dist = 5, match_no_auth = TRUE, match_canon = TRUE) %>%
    as_tibble(),
  match_results_classified_round_3 = tt_classify_result(match_results_raw_round_3),
  match_results_resolved_round_3 = tt_resolve_synonyms(match_results_classified_round_3, world_ferns_data),
  # Combine name resolution results
  match_results_resolved_all =
    combined_match_results(
      ncbi_names_query = ncbi_names_query, 
      match_results_resolved_round_1, match_results_resolved_round_2, match_results_resolved_round_3),
  # Map NCBI names to accepted names
  ncbi_accepted_names_map = make_ncbi_accepted_names_map(match_results_resolved_all),
  
  # Remove rogue sequences ----
  # Combine sanger sequences and metadata, filter to resolved names
  sanger_seqs_combined_filtered = combine_and_filter_sanger(raw_meta, raw_fasta, ncbi_accepted_names_map, ppgi_taxonomy),
  # Conduct all-by-all blast
  all_by_all_blast = blast_rogues(sanger_seqs_combined_filtered),
  # Identify rogues (sequences matching wrong family)
  # FIXME: many of these are due to bad taxonomy. inspect results and modify WF taxonomy as needed
  sanger_seqs_rogues = detect_rogues(
    metadata_with_seqs = sanger_seqs_combined_filtered,
    blast_results = all_by_all_blast,
    ppgi = ppgi_taxonomy,
    id = gene),
  sanger_seqs_rogues_removed = anti_join(
    sanger_seqs_combined_filtered,
    sanger_seqs_rogues,
    by = c("accession", "gene")
  ),

  # Select final Sanger sequences ----
  # Select one specimen per species, prioritizing in order
  # - 1: specimens with rbcL + any other gene
  # - 2: specimens with rbcL
  # - 3: specimens with longest combined non-rbcL genes
  sanger_accessions_selection = select_genbank_genes(sanger_seqs_rogues_removed)
)
