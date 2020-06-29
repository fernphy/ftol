# Define variables used to generate plan
#
# - Target plastid fern genes to download individually 
# (not as part of a plastome)
target_genes <- c("rbcL", "atpA", "atpB", "rps4")
# - Minimum lengths for each gene (in case this needs to be adjusted per gene)
min_lengths <- c(400, 400, 400, 400)
# - Most recent date cutoff for sampling genes
date_cutoff <- "2020/06/28"

plan <- drake_plan(
  
  # Load data ----
  
  # Use `target("path", format = "file")` to track file contents
  
  # Data for resolving taxonomic names:
  # Catalog of Life database subsetted to tracheophytes.
  col_plants_path = target("data_raw/archive-kingdom-plantae-phylum-tracheophyta-bl3/taxa.txt", format = "file"),
  col_plants = pferns::load_col_plants(col_plants_path),
  
  # Modified PPGI taxonomy
  # with new genera and slightly different treatments following World Ferns list
  ppgi_taxonomy_path = target("data_raw/ppgi_taxonomy_mod.csv", format = "file"),
  ppgi_taxonomy = read_csv(ppgi_taxonomy_path),
  
  # List of plastid coding genes from Wei et al 2017
  target_sanger_genes_path = target("data_raw/wei_2017_coding_genes.txt", format = "file"),
  target_sanger_genes = read_lines(target_sanger_genes_path),
  
  # Outgroup plastome accessions
  plastome_outgroups_path = target("data_raw/plastome_outgroups.csv", format = "file"),
  plastome_outgroups = read_csv(plastome_outgroups_path),
  
  # Calibration dates after Testo and Sundue 2016
  plastome_calibration_dates_path = target("data_raw/testo_sundue_2016_calibrations.csv", format = "file"),
  plastome_calibration_dates = load_calibration_dates(plastome_calibration_dates_path),
  
  # Manually selected synonyms for resolving names of plastid genes
  genbank_names_with_mult_syns_select_path = target("data_raw/genbank_names_with_mult_syns_select.csv", format = "file"),
  genbank_names_with_mult_syns_select = read_csv(genbank_names_with_mult_syns_select_path),
  
  # Download individual plastid sequences from GenBank----
  
  # Download fern plastid gene fasta files
  # (target genes defined in plastid_make.R)
  raw_fasta = target(
    fetch_fern_gene(gene, end_date = date_cutoff),
    transform = map(gene = !!target_genes)
  ),
  
  # Download fern plastid gene metadata
  raw_meta = target(
    fetch_fern_metadata(gene, end_date = date_cutoff),
    transform = map(gene = !!target_genes)
  ),
  
  # Combine GenBank sequences with metadata 
  # (so fasta sequence is a list-column in genbank_seqs_combined_raw)
  genbank_seqs_combined_raw = target(
    join_genbank_fasta_with_meta(
      seqs = raw_fasta, 
      metadata = raw_meta, 
      id = gene),
    transform = map(
      raw_fasta, 
      raw_meta, 
      gene = !!target_genes, 
      .id = gene)
  ),
  
  # Filter by minimum length per gene and genus in PPGI
  genbank_seqs_combined_filtered = target(
    filter_genbank_seqs(
      metadata_with_seqs = genbank_seqs_combined_raw, 
      min_len = min_len, 
      ppgi = ppgi_taxonomy,
      id = gene),
    transform = map(
      genbank_seqs_combined_raw,
      min_len = !!min_lengths,
      gene = !!target_genes, 
      .id = gene)
  ),
  
  # Filter out rogue sequences
  # - run all-by-all BLAST to identify sequences that match the wrong family
  all_by_all_blast = target(
    blast_rogues(
      metadata_with_seqs = genbank_seqs_combined_filtered,
      id = gene),
    transform = map(
      genbank_seqs_combined_filtered,
      gene = !!target_genes, 
      .id = gene)
  ),
  
  # - identify rogue sequences
  genbank_seqs_rogues = target(
    detect_rogues(
      metadata_with_seqs = genbank_seqs_combined_filtered, 
      blast_results = all_by_all_blast,
      ppgi = ppgi_taxonomy,
      id = gene),
    transform = map(
      genbank_seqs_combined_filtered, 
      all_by_all_blast,
      gene = !!target_genes,
      .id = gene)
  ),
  
  # - remove rogue sequences
  genbank_seqs_rogues_removed = target(
    anti_join_tracked(
      id = gene,
      genbank_seqs_combined_filtered,
      genbank_seqs_rogues,
      by = "accession"),
    transform = map(
      genbank_seqs_combined_filtered,
      genbank_seqs_rogues,
      gene = !!target_genes,
      .id = gene)
  ),
  
  # Combine cleaned genes into single dataframe
  genbank_seqs_rogues_removed_all_genes = target(
    bind_rows(genbank_seqs_rogues_removed, .id = "gene"),
    transform = combine(genbank_seqs_rogues_removed, .id = gene),
  ),

  # Resolve names:
  # - first run automatic name resolution
  genbank_seqs_names_resolved_auto = resolve_genbank_names_auto(
    genbank_seqs_rogues_removed_all_genes,
    col_plants
  ),

  # (manually select synonyms, save as genbank_names_with_mult_syns_select.csv)

  # - then determine final names to use, filtering out names that didn't match.
  genbank_seqs_names_resolved = resolve_genbank_names_final(
    names_resolved_auto = genbank_seqs_names_resolved_auto$names_resolved_auto,
    names_resolved_to_other_sources = genbank_seqs_names_resolved_auto$names_resolved_to_other_sources,
    name_resolution_syns_to_use = genbank_names_with_mult_syns_select,
    combined_metadata = genbank_seqs_rogues_removed_all_genes
  ),

  # Select final GenBank accessions.
  # Select one specimen per species, prioritizing in order
  # - 1: specimens with rbcL + any other gene
  # - 2: specimens with rbcL
  # - 3: specimens with longest combined non-rbcL genes
  genbank_accessions_selection = select_genbank_genes(
    genbank_seqs_names_resolved, target_genes),
  
  # Download core set of plastid genes from plastomes ----
  # ca. 100 species by 60 genes

  # Download plastome metadata (accessions and species)
  plastome_metadata = download_plastome_metadata(
    end_date = date_cutoff,
    outgroups = plastome_outgroups),

  # Resolve species names in metadata using CoL plants as the taxonomic standard.
  plastome_metadata_renamed = resolve_pterido_plastome_names(
    plastome_metadata,
    col_plants),

  # Extract list of accessions for looping
  plastome_accessions = plastome_metadata_renamed$accession,

  # Download the coding genes for each plastome one at a time.
  plastid_seqs = target(
    gbfetch::fetch_gene_from_genome(
      target_sanger_genes,
      plastome_accessions),
    dynamic = map(plastome_accessions)
  ),

  # Clean up the downloaded sequences: remove NULL values, duplicated genes, etc.
  cleaned_plastid_seqs = target(
    clean_plastid_seqs(plastid_seqs),
    dynamic = map(plastid_seqs)
  ),

  # Combine results of cleaning downloaded sequences
  # and set names by accession.
  # Use readd() to convert the object back to a static object from dynamic.
  cleaned_plastid_seqs_list = readd(cleaned_plastid_seqs, cache = plastid_cache) %>%
    purrr::set_names(plastome_accessions),

  # Select final accessions / genes
  # - best representative accession per species
  # - only include genes and accessions with > 50% occupancy.
  plastid_selection = select_plastid_seqs(
    cleaned_plastid_seqs_list, plastome_metadata_renamed, "species"),

  # Reformat as list of unaligned genes.
  plastid_genes_unaligned = extract_seqs_by_gene(cleaned_plastid_seqs_list, plastid_selection),

  # Combine genes from GenBank with genes from plastomes ----

  # Combine raw GenBank fasta sequences.
  raw_fasta_all_genes = target(
    bind_rows(genbank_seqs_combined_raw, .id = "gene"),
    transform = combine(genbank_seqs_combined_raw),
  ),

  # Combine genes from GenBank with genes from plastomes (still unaligned).
  plastid_genes_unaligned_combined = combine_genbank_with_plastome(
    raw_fasta_all_genes,
    target_genes,
    genbank_accessions_selection,
    plastome_metadata_renamed,
    plastid_genes_unaligned
  ),

  # Align each gene.
  plastid_genes_aligned = purrr::map(
    plastid_genes_unaligned_combined,
    ~ips::mafft(
      x = .,
      options = "--adjustdirection",
      exec = "/home/nittaj/miniconda3/envs/pacferns/bin/mafft")),

  # Trim alignments.
  plastid_genes_aligned_trimmed = purrr::map(
    plastid_genes_aligned,
    trimal_auto),

  # Rename sequences in each gene as species
  # (using underscores, not spaces).
  # - first combine accession + resolved name for GenBank genes and
  # plastome genes
  resolved_names_all = bind_rows(
    select(genbank_seqs_names_resolved, accession, species),
    select(plastome_metadata_renamed, accession, species)
  ) %>% unique,

  # - then use the combined, resolved names to rename each alignment
  # by species instead of accession
  plastid_genes_aligned_trimmed_renamed = purrr::map(
    plastid_genes_aligned_trimmed, rename_alignment,
    name_metadata = resolved_names_all),

  # Concatenate alignments by species name.
  plastome_alignment = concatenate_genes(
    plastid_genes_aligned_trimmed_renamed),

  # Generate tree.
  plastome_tree = jntools::iqtree(
    plastome_alignment,
    m = "GTR+I+G", bb = 1000, nt = "AUTO",
    redo = FALSE, echo = TRUE, wd = here::here("iqtree")),

  # Dating analysis with treepl ----

  # Root tree on bryophytes
  plastome_tree_rooted = ape::root(
    plastome_tree,
    c("Anthoceros_angustus", "Marchantia_polymorpha", "Physcomitrella_patens")),

  # Run initial treepl search to identify smoothing parameter
  treepl_cv_results = run_treepl_cv(
    phy = plastome_tree_rooted,
    alignment = plastome_alignment,
    calibration_dates = plastome_calibration_dates,
    cvstart = "1000",
    cvstop = "0.000001",
    plsimaniter = "200000", # preliminary output suggested > 100000
    seed = 7167,
    thorough = TRUE,
    wd = here::here("treepl"),
    nthreads = 1,
    echo = TRUE
  ),

  # Run priming analysis to determine optimal states for other parameters
  treepl_priming_results = run_treepl_prime(
    phy = plastome_tree_rooted,
    alignment = plastome_alignment,
    calibration_dates = plastome_calibration_dates,
    cv_results = treepl_cv_results,
    plsimaniter = "200000", # preliminary output suggested > 100000
    seed = 7167,
    thorough = TRUE,
    wd = here::here("treepl"),
    nthreads = 1,
    echo = TRUE
  ),

  # Run treePL dating analysis
  treepl_dating_results = run_treepl(
    phy = plastome_tree_rooted,
    alignment = plastome_alignment,
    calibration_dates = plastome_calibration_dates,
    cv_results = treepl_cv_results,
    priming_results = treepl_priming_results,
    plsimaniter = "200000", # preliminary output suggested > 100000
    seed = 7167,
    thorough = TRUE,
    wd = here::here("treepl"),
    nthreads = 7,
    echo = TRUE
  )

)
