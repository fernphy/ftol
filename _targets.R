library(targets)
library(tarchetypes)
source("R/packages.R")
source("R/functions.R")

# Specify path to folder with raw data
data_raw <- "/data_raw"

# Set parallel back-end
plan(callr)

# Use targets workspaces for debugging
tar_option_set(workspace_on_error = TRUE)
# Return tibble from taxastand functions by default
options(ts_tbl_out = TRUE)

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
  # List of coding genes to extract from plastomes
  # (based on genes of Wei et al 2017, then trimmed to non-duplicated genes)
  tar_file(target_plastome_genes_path, path(data_raw, "target_coding_genes.txt")),
  target_plastome_genes = read_lines(target_plastome_genes_path),
  # Outgroup plastome accessions
  tar_file(plastome_outgroups_path, path(data_raw, "plastome_outgroups.csv")),
  plastome_outgroups = read_csv(plastome_outgroups_path),
  # Calibration dates after Testo and Sundue 2016
  tar_file(plastome_calibration_dates_path, path(data_raw, "testo_sundue_2016_calibrations.csv")),
  plastid_calibration_dates = load_calibration_dates(plastome_calibration_dates_path),
  # Manually selected synonyms for resolving names of plastid genes
  tar_file(sanger_names_with_mult_syns_select_path, path(data_raw, "genbank_names_with_mult_syns_select.csv")),
  sanger_names_with_mult_syns_select = read_csv(sanger_names_with_mult_syns_select_path),
  # Manually curated list of raw fasta accessions to exclude from analysis
  tar_file(raw_fasta_exclude_list_path, path(data_raw, "raw_fasta_exclude_list.csv")),
  raw_fasta_exclude_list = read_csv(raw_fasta_exclude_list_path),

  # Prep for assembling Sanger plastid regions ----
  # Define variables used in plan:
  # - Target plastic loci (coding genes and spacers)
  target_loci = c("atpA", "atpB", "matK", "rbcL", "rps4", "trnL-trnF", "rps4-trnS"),
  # - Most recent date cutoff for sampling genes
  date_cutoff = "2021/10/26",

  # Assemble reference Sanger sequences ----
  # Download raw gene sequences for each target locus
  tar_target(
    fern_ref_seqs_raw,
    fetch_fern_ref_seqs(
      target = target_loci, end_date = date_cutoff, 
      # Only include intron for trnL-trnF
      req_intron = str_detect(target_loci, "trnL-trnF"), 
      workers = 48),
    pattern = map(target_loci),
    deployment = "main"
  ),
  # Filter to one longest seq per genus per target
  fern_ref_seqs = filter_ref_seqs(fern_ref_seqs_raw),
  # Align sequences
  tar_target(
    fern_ref_seqs_aligned,
    align_ref_seqs(fern_ref_seqs, target_loci),
    pattern = map(target_loci)
  ),
  # Trim sequences (lightly)
  # - dnabin format for IQTREE
  tar_target(
    fern_ref_seqs_trimmed_dnabin,
    trimal(
      fern_ref_seqs_aligned,
      c("-gt", "0.05", "-terminalonly"),
      return_seqtbl = FALSE),
    pattern = map(fern_ref_seqs_aligned)
  ),
  # - seqtbl format for downstream steps
  tar_target(
    fern_ref_seqs_trimmed,
    dnabin_to_seqtbl(fern_ref_seqs_trimmed_dnabin),
    pattern = map(fern_ref_seqs_trimmed_dnabin)
  ),
  # Convert to format for extract_from_ref_blast()
  fern_ref_seqs_trimmed_clean = separate(
    fern_ref_seqs_trimmed,
    accession,
    c("accession", "species", "target"),
    sep = "__"),
  # Build trees (to check quality of reference sequences)
  tar_target(
    fern_ref_seqs_tree,
    jntools::iqtree(
      fern_ref_seqs_trimmed_dnabin,
      m = "GTR+I+G", bb = 1000, nt = "AUTO",
      redo = TRUE, echo = TRUE, wd = tempdir()
      ),
    pattern = map(fern_ref_seqs_trimmed_dnabin)
  ),
  # Write out alignments for inspection
  tar_file(
    fern_ref_seqs_trimmed_out,
    write_fasta_tar(
      fern_ref_seqs_trimmed_dnabin,
      paste0("results/ref_seqs/", target_loci, "_ref_aln.fasta")
    ),
    pattern = map(fern_ref_seqs_trimmed_dnabin, target_loci)
  ),
  # Write out trees for inspection
  tar_file(
    fern_ref_seqs_tree_out,
    write_tree_tar(
      fern_ref_seqs_tree,
      paste0("results/ref_seqs/", target_loci, "_ref_phy.tree")
    ),
    pattern = map(fern_ref_seqs_tree, target_loci)
  ),

  # Download individual plastid sequences from GenBank (Sanger sequences) ----
  
  # Download raw fasta files
  tar_target(
    fern_sanger_seqs_raw,
    fetch_fern_sanger_seqs(target_loci, end_date = date_cutoff, accs_exclude_list = NULL),
    pattern = map(target_loci),
    deployment = "main" # don't run in parallel, or will get HTTP status 429 errors
  ),
  # Extract target regions
  flavors = c("blastn", "dc-megablast"),
  tar_target(
    fern_sanger_extract_res,
    extract_from_ref_blast(
      query_seqtbl = fern_sanger_seqs_raw,
      ref_seqtbl = fern_ref_seqs_trimmed_clean,
      target = target_loci,
      blast_flavor = flavors,
      other_args = c("-m", "span", "--threads", "4")
    ),
    pattern = cross(target_loci, flavors)
  ),
  # after comparing results between blastn and dc-megablast, dc-megablast seems to work better
  raw_fasta = clean_extract_res(fern_sanger_extract_res, "dc-megablast"),
  # Fetch metadata
  tar_target(
    raw_meta_all,
    fetch_fern_metadata(target_loci, end_date = date_cutoff),
    pattern = map(target_loci),
    deployment = "main"
  ),
  raw_meta = unique(raw_meta_all),
  
  # Taxonomic name resolution ----
  # Download species names from NCBI
  ncbi_names_raw = raw_meta %>% pull(taxid) %>% unique %>% fetch_taxonomy,
  # Clean NCBI species names
  ncbi_names_full = clean_ncbi_names(ncbi_names_raw),
  # Exclude invalid names (hybrids, taxa not identified to species level)
  ncbi_names_query = exclude_invalid_ncbi_names(ncbi_names_full),
  # Parse reference names
  wf_ref_names = ts_parse_names(unique(world_ferns_data$scientificName)),
  # Resolve names, round 1: NCBI accepted scientific names
  ncbi_names_query_round_1 = select_ncbi_names_round_1(ncbi_names_query),
  # - match names to reference
  match_results_raw_round_1 = ts_match_names(
    query = ncbi_names_query_round_1$scientific_name, 
    reference = wf_ref_names,
    max_dist = 5, match_no_auth = TRUE, match_canon = TRUE),
  # - resolve synonyms
  match_results_resolved_round_1 = ts_resolve_names(match_results_raw_round_1, world_ferns_data),
  # Resolve names, round 2: NCBI synonym scientific names
  ncbi_names_query_round_2 = select_ncbi_names_round_2(match_results_resolved_round_1, ncbi_names_query),
  match_results_raw_round_2 = ts_match_names(
    query = ncbi_names_query_round_2$scientific_name, 
    reference = wf_ref_names,
    max_dist = 5, match_no_auth = TRUE, match_canon = TRUE),
  match_results_resolved_round_2 = ts_resolve_names(match_results_raw_round_2, world_ferns_data),
  # Resolve names, round 3: NCBI species without author
  ncbi_names_query_round_3 = select_ncbi_names_round_3(match_results_resolved_round_1, match_results_resolved_round_2, ncbi_names_query),
  match_results_raw_round_3 = ts_match_names(
    query = ncbi_names_query_round_3$species, 
    reference = wf_ref_names,
    max_dist = 5, match_no_auth = TRUE, match_canon = TRUE),
  match_results_resolved_round_3 = ts_resolve_names(match_results_raw_round_3, world_ferns_data),
  # Combine name resolution results
  match_results_resolved_all =
    combined_match_results(
      ncbi_names_query = ncbi_names_query, 
      match_results_resolved_round_1, match_results_resolved_round_2, match_results_resolved_round_3),
  # Map NCBI names to accepted names
  ncbi_accepted_names_map = make_ncbi_accepted_names_map(match_results_resolved_all),
  
  # Remove rogue sequences ----
  # Combine sanger sequences and metadata, filter to resolved names
  # - set minimum lengths (bp) for filtering genes and spacers
  min_gene_len = 200,
  min_spacer_len = 20,
  sanger_seqs_combined_filtered = combine_and_filter_sanger(
    raw_meta, raw_fasta, ncbi_accepted_names_map, 
    ppgi_taxonomy, min_gene_len, min_spacer_len),
  # Make BLAST database including all fern sequences
  tar_file(
    sanger_blast_db,
    make_fern_blast_db(
      seqtbl = sanger_seqs_combined_filtered, 
      blast_db_dir = "intermediates/blast_sanger", 
      out_name = "ferns_sanger")
  ),
  # Group query sequences for parallel computing
  tar_group_count(
    sanger_blast_query,
    dplyr::select(sanger_seqs_combined_filtered, seq, otu),
    count = 30), # number of jobs to run in parallel
  # Conduct all-by-all blast in parallel
  tar_target(
    all_by_all_blast,
    blast_seqtbl(
      seqtbl = sanger_blast_query,
      blastdb_files = sanger_blast_db
    ),
    pattern = map(sanger_blast_query)
  ),
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
  sanger_accessions_selection = select_genbank_genes(sanger_seqs_rogues_removed),

  # Download core set of plastid genes from plastomes ----
  # Download plastome metadata (accessions and species)
  plastome_metadata_raw = download_plastome_metadata(
    end_date = date_cutoff,
    outgroups = plastome_outgroups),
  # Resolve species names in plastome metadata
  # (will drop accession if name could not be resolved)
  # FIXME: add Lepisorus hederaceus (Christ) R.Wei & X.C.Zhang
  # and Elaphoglossum marginatum var. marginatum to World Ferns taxonomy
  plastome_metadata_raw_renamed = resolve_pterido_plastome_names(
    plastome_metadata_raw, plastome_outgroups, wf_ref_names, world_ferns_data
  ),
  # Download plastome sequences
  # don't run in parallel, or will get HTTP status 429 errors
  target_plastome_accessions = unique(plastome_metadata_raw_renamed$accession),
  tar_target(
    plastome_genes_raw,
    fetch_fern_genes_from_plastome(
      genes = target_plastome_genes, 
      accession = target_plastome_accessions),
    pattern = map(target_plastome_accessions),
    deployment = "main"),
  # Combine plastome metadata and sequences, filter to best accession per taxon
  plastome_seqs_combined_filtered = select_plastome_seqs(
    plastome_genes_raw, plastome_metadata_raw_renamed),

  # Combine and align Sanger and plastome sequences ----
  # Combine Sanger and plastome sequences into single dataframe, group by gene
  tar_group_by(
    plastid_genes_unaligned_combined,
    combine_sanger_plastome(
      sanger_accessions_selection,
      sanger_seqs_combined_filtered,
      plastome_seqs_combined_filtered),
    gene),
  # Align sequences by gene
  tar_target(
    plastid_genes_aligned,
    align_seqs_tbl(plastid_genes_unaligned_combined),
    pattern = map(plastid_genes_unaligned_combined)
  ),
  # Trim alignments
  tar_target(
    plastid_genes_aligned_trimmed,
    mutate(
      plastid_genes_aligned,
      seq = purrr::map(seq, ~trimal(., other_args = "-automated1"))
    ),
    pattern = map(plastid_genes_aligned)
  ),
  # Concatenate alignments
  plastid_alignment = do.call(
    ape::cbind.DNAbin, 
    c(plastid_genes_aligned_trimmed$seq, fill.with.gaps = TRUE)
  ),

  # Phylogenetic analysis
  # Generate tree: single concatenated analysis.
  plastid_tree = jntools::iqtree(
    plastid_alignment,
    m = "GTR+I+G", bb = 1000, nt = "AUTO",
    redo = TRUE, echo = TRUE, wd = here::here("intermediates/iqtree")
  )

)
