# Define variables used to generate plan
#
# - Target plastid fern genes to download individually
# (not as part of a plastome)
target_genes <- c("rbcL", "atpA", "atpB", "rps4")
# - Minimum lengths for each gene (in case this needs to be adjusted per gene)
min_lengths <- c(400, 400, 400, 400)
# - Most recent date cutoff for sampling genes
date_cutoff <- "2020/06/30"

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
  plastid_calibration_dates = load_calibration_dates(plastome_calibration_dates_path),

  # Manually selected synonyms for resolving names of plastid genes
  sanger_names_with_mult_syns_select_path = target("data_raw/genbank_names_with_mult_syns_select.csv", format = "file"),
  sanger_names_with_mult_syns_select = read_csv(sanger_names_with_mult_syns_select_path),

  # GoFlag (short reads from seq-capture and genome skimming) metadata
  goflag_meta_path = target("data_raw/goflag/Pilot_Ferns_TargetCapture_Skimming.txt", format = "file"),
  goflag_meta = read_tsv(goflag_meta_path) %>%
    select(taxon = Taxon, targeted_capture_id = `Targeted Capture ID`, genome_skimming_id  = `Genome Skimming ID`) %>%
    filter(!is.na(taxon)),

  # Download individual plastid sequences from GenBank (Sanger sequences) ----

  # Download fern plastid gene fasta files
  # (target genes defined in plastid_make.R)
  raw_fasta = target(
    fetch_fern_gene(gene, end_date = date_cutoff),
    transform = map(gene = !!target_genes),
    hpc = FALSE # don't run in parallel, or will get HTTP status 429 errors
  ),

  # Download fern plastid gene metadata
  raw_meta = target(
    fetch_fern_metadata(gene, end_date = date_cutoff),
    transform = map(gene = !!target_genes),
    hpc = FALSE
  ),

  # Standardize names and filter Sanger sequences ----

  # Combine GenBank sequences with metadata
  # (so fasta sequence is a list-column in sanger_seqs_combined_raw)
  sanger_seqs_combined_raw = target(
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
  sanger_seqs_combined_filtered = target(
    filter_sanger_seqs(
      metadata_with_seqs = sanger_seqs_combined_raw,
      min_len = min_len,
      ppgi = ppgi_taxonomy,
      id = gene),
    transform = map(
      sanger_seqs_combined_raw,
      min_len = !!min_lengths,
      gene = !!target_genes,
      .id = gene)
  ),

  # Filter out rogue sequences
  # - run all-by-all BLAST to identify sequences that match the wrong family
  all_by_all_blast = target(
    blast_rogues(
      metadata_with_seqs = sanger_seqs_combined_filtered,
      id = gene),
    transform = map(
      sanger_seqs_combined_filtered,
      gene = !!target_genes,
      .id = gene)
  ),

  # - identify rogue sequences
  sanger_seqs_rogues = target(
    detect_rogues(
      metadata_with_seqs = sanger_seqs_combined_filtered,
      blast_results = all_by_all_blast,
      ppgi = ppgi_taxonomy,
      id = gene),
    transform = map(
      sanger_seqs_combined_filtered,
      all_by_all_blast,
      gene = !!target_genes,
      .id = gene)
  ),

  # - remove rogue sequences
  sanger_seqs_rogues_removed = target(
    anti_join_tracked(
      id = gene,
      sanger_seqs_combined_filtered,
      sanger_seqs_rogues,
      by = "accession"),
    transform = map(
      sanger_seqs_combined_filtered,
      sanger_seqs_rogues,
      gene = !!target_genes,
      .id = gene)
  ),

  # Combine cleaned genes into single dataframe
  sanger_seqs_rogues_removed_all_genes_raw = target(
    bind_rows(sanger_seqs_rogues_removed, .id = "gene"),
    transform = combine(sanger_seqs_rogues_removed, .id = gene),
  ),

  # - convert genes from number codes to gene names
  sanger_seqs_rogues_removed_all_genes = rename_genes(
    sanger_seqs_rogues_removed_all_genes_raw,
    target_genes
  ),

  # Resolve names:
  # - first run automatic name resolution
  sanger_seqs_names_resolved_auto = resolve_sanger_names_auto(
    sanger_seqs_rogues_removed_all_genes,
    col_plants
  ),

  # (manually select synonyms, save as sanger_names_with_mult_syns_select.csv)

  # - then determine final names to use, filtering out names that didn't match.
  sanger_seqs_names_resolved = resolve_sanger_names_final(
    names_resolved_auto = sanger_seqs_names_resolved_auto$names_resolved_auto,
    names_resolved_to_other_sources = sanger_seqs_names_resolved_auto$names_resolved_to_other_sources,
    name_resolution_syns_to_use = sanger_names_with_mult_syns_select,
    combined_metadata = sanger_seqs_rogues_removed_all_genes
  ),

  # Select final Sanger sequences (GenBank accessions).
  #
  # Select one specimen per species, prioritizing in order
  # - 1: specimens with rbcL + any other gene
  # - 2: specimens with rbcL
  # - 3: specimens with longest combined non-rbcL genes
  sanger_accessions_selection = select_genbank_genes(
      genbank_seqs_tibble = sanger_seqs_names_resolved,
      n_seqs_per_sp = "single"
  ),

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

  # Extract the coding genes for each plastome one at a time.
  plastome_seqs = target(
    fetch_fern_genes_from_plastome(
      target_sanger_genes,
      plastome_accessions),
    dynamic = map(plastome_accessions)
  ),

  # Combine plastome sequences and set names by accession.
  # Use readd() to convert the object back to a static object from dynamic.
  plastome_seqs_list = readd(plastome_seqs, cache = plastid_cache),

  # Select final accessions / genes, only include genes and accessions with > 50% occupancy
  plastome_selection = select_plastome_seqs(
    plastome_seqs_list,
    plastome_metadata_renamed,
    filter_by = "species"
    ),

  # Reformat as list of unaligned genes.
  plastome_genes_unaligned = extract_seqs_by_gene(plastome_seqs_list, plastome_selection),

  # Skim plastid genes from short-read data ----

  # Part 1: Trim raw fastq files with fastp

  # Combine GoFlag target capture and skimming samples into single vector
  seq_cap_samples = c(goflag_meta$targeted_capture_id, goflag_meta$genome_skimming_id),

  # Run fastp on each sample.
  # Writes trimmed fastq files to intermediates/fastp/
  # and returns a dataframe with trimming stats.
  trim_results_each = target(
    fastp(seq_cap_samples),
    dynamic = map(seq_cap_samples)
  ),

  # Combine results into single summary dataframe
  trim_results_summary = target(
    bind_rows(trim_results_each),
    dynamic = group(trim_results_each)
  ),

  # Part 2: Format "targets" (genes to extract from short read data)

  # Assemble set of coding genes from GenBank plastome data
  # using map(accessions, ~fetch_genes_from_plastome(., target_sanger_genes))
  # Run this command locally as juno seems to have a hard time keeping a connection
  hybpiper_plastid_targets_path = target("temp/plastid_targets.RDS", format = "file"),
  hybpiper_plastid_targets = readRDS(hybpiper_plastid_targets_path),

  # Collapse amino acid targets into single list, write it out
  hybpiper_aa_targets_out = ape::write.FASTA(
    transpose(hybpiper_plastid_targets)[["aa"]] %>% do.call(c, .),
    file_out("intermediates/hybpiper/plastid_aa_targets.fasta")
  ),

  # Part 3: Run HybPiper to sort reads into genes

  # Get vector of trimmed reads for hybpiper
  forward_reads = get_reads(
    data_dir = here::here("intermediates/fastp/"),
    pattern = "R1.fastq",
    depends = trim_results_summary
  ),

  reverse_reads = get_reads(
    data_dir = here::here("intermediates/fastp/"),
    pattern = "R2.fastq",
    depends = trim_results_summary
  ),

  # Make list of paired reads for HybPiper
  paired_reads_list = make_paired_reads_list (
    forward_reads,
    reverse_reads
  ),

  # Map HybPiper 'reads_first' over the reads by sample.
  # Note (with blastx) it takes ca. 28 min per sample with 1 CPU, 18 min with 2 CPU,
  # 13 min with 6 CPU, and 7 min with 10 CPU
  # So for a large number of samples (eg 40) it is better to run 40 CPU in
  # parallel on 40 samples with 1 CPU each
  hybpiper_results_each = target(
    reads_first(
      wd = here::here("intermediates/hybpiper"),
      echo = FALSE,
      # Use amino-acids baitfile
      baitfile = file_in("intermediates/hybpiper/plastid_aa_targets.fasta"),
      # When paired_reads_list gets split up by dynamic mapping, it is split into lists.
      # The character vector we want for `readfiles` is the first element of each list
      readfiles = paired_reads_list[[1]],
      prefix = paired_reads_list[[1]][[1]] %>% fs::path_file() %>% str_remove_all("_R.\\.fastq"),
      cpu = 1,
      # use blastx, not BWA
      bwa = FALSE,
      # don't run exonerate or try to assemble genes
      other_args = c("--no-exonerate", "--no-assemble")),
    dynamic = map(paired_reads_list)
  ),

  # Combine results
  hybpiper_results = target(
    c(hybpiper_results_each),
    dynamic = group(hybpiper_results_each)
  ),

  # Part 4: align reads to reference, extract consensus

  # Align read fragments to reference
  # - loop over the list of samples
  each_extracted_reads_consensus = target(
    get_hybpip_consensus(
      sample = seq_cap_samples,
      plastid_targets = hybpiper_plastid_targets,
      depends = hybpiper_results),
    dynamic = map(seq_cap_samples)
  ),

  # - combine results
  extracted_reads_consensus = target(
    bind_rows(each_extracted_reads_consensus),
    dynamic = group(each_extracted_reads_consensus)
  ),

  # Combine datasets ----

  # Combine raw GenBank (Sanger) fasta sequences
  raw_fasta_all_genes = target(
    bind_rows(sanger_seqs_combined_raw, .id = "gene"),
    transform = combine(sanger_seqs_combined_raw),
  ),

  # Rename genes as gene names (not numbers)
  raw_fasta_all_genes_renamed = rename_genes(
    raw_fasta_all_genes,
    target_genes
  ),

  # Filter out species already in plastomes from Sanger sequences
  sanger_accessions_selection_filtered =
    filter_out_plastome_species(
      plastome_genes_unaligned = plastome_genes_unaligned,
      plastome_metadata_renamed = plastome_metadata_renamed,
      sanger_accessions_selection = sanger_accessions_selection,
      filter = TRUE
    ),

  # Combine genes from GenBank with genes from plastomes (still unaligned).
  plastid_genes_unaligned_combined =
    combine_sanger_with_plastome(
      raw_fasta_all_genes = raw_fasta_all_genes_renamed,
      sanger_accessions_selection_filtered,
      plastome_genes_unaligned
  ),

  # Align each gene.
  plastid_genes_aligned = purrr::map(
      plastid_genes_unaligned_combined,
      ~ips::mafft(
        x = .,
        options = "--adjustdirection",
        exec = "/usr/bin/mafft")
  ),

  # Trim alignments.
  plastid_genes_aligned_trimmed = purrr::map(
      plastid_genes_aligned,
      trimal_auto
  ),

  # Rename sequences in each gene as species
  # (using underscores, not spaces).
  # - first combine accession + resolved name for GenBank genes and
  # plastome genes
  resolved_names_all = bind_rows(
    select(sanger_seqs_names_resolved, accession, species),
    select(plastome_metadata_renamed, accession, species)
  ) %>% unique,

  # - then use the combined, resolved names to rename each alignment
  # by species instead of accession
  plastid_genes_aligned_trimmed_renamed = rename_alignment_list(
    alignment_list = plastid_genes_aligned_trimmed,
    name_metadata = resolved_names_all
  ),

  # Concatenate alignments by species name.
  plastid_alignment = concatenate_genes(plastid_genes_aligned_trimmed_renamed),

  # Generate tree.
  plastid_tree = jntools::iqtree(
     plastid_alignment,
     m = "GTR+I+G", bb = 1000, nt = "AUTO",
     redo = FALSE, echo = TRUE, wd = here::here("intermediates/iqtree")
  ),

  # Dating analysis with treepl ----

  # Root tree on bryophytes
  plastid_tree_rooted = ape::root(
      plastid_tree,
      c("Anthoceros_angustus", "Marchantia_polymorpha", "Physcomitrium_patens")),

  # Run initial treepl search to identify smoothing parameter
  treepl_cv_results = run_treepl_cv(
      phy = plastid_tree_rooted,
      alignment = plastid_alignment,
      calibration_dates = plastid_calibration_dates,
      cvstart = "1000",
      cvstop = "0.000001",
      plsimaniter = "200000", # preliminary output suggested > 100000
      seed = 7167,
      thorough = TRUE,
      wd = here::here("intermediates/treepl"),
      nthreads = 1,
      echo = TRUE
    ),

  # Run priming analysis to determine optimal states for other parameters
  treepl_priming_results = run_treepl_prime(
      phy = plastid_tree_rooted,
      alignment = plastid_alignment,
      calibration_dates = plastid_calibration_dates,
      cv_results = treepl_cv_results,
      plsimaniter = "200000", # preliminary output suggested > 100000
      seed = 7167,
      thorough = TRUE,
      wd = here::here("intermediates/treepl"),
      nthreads = 1,
      echo = TRUE
    ),

  # Run treePL dating analysis
  plastid_tree_dated = run_treepl(
      phy = plastid_tree_rooted,
      alignment = plastid_alignment,
      calibration_dates = plastid_calibration_dates,
      cv_results = treepl_cv_results,
      priming_results = treepl_priming_results,
      plsimaniter = "200000", # preliminary output suggested > 100000
      seed = 7167,
      thorough = TRUE,
      wd = here::here("intermediates/treepl"),
      nthreads = 7,
      echo = TRUE
    ),

  # Output trees and alignments to results
  
  # - write outplastid tree (not dated)
  plastid_tree_out = ape::write.tree(
    plastid_tree_rooted, 
    file_out("results/releases/ftol_plastid.tre")),

  # - write outplastid tree (dated)
  plastid_tree_dated_out = ape::write.tree(
    plastid_tree_dated, 
    file_out("results/releases/ftol_plastid_dated.tre")),

  # - write out plastid concatenated alignment
  plastid_alignment_concat_out = ape::write.FASTA(
    plastid_alignment, 
    file_out("results/releases/ftol_plastid_concat.fasta")),

  # - make table of GenBank accession numbers 
  plastid_acc_data = make_acc_ref_table(
    plastid_genes_aligned_trimmed = plastid_genes_aligned_trimmed,
    sanger_seqs_names_resolved = sanger_seqs_names_resolved,
    plastome_metadata_renamed = plastome_metadata_renamed,
    target_genes = target_genes),
  
  # - write out table of GenBank accession numbers 
  plastid_acc_data_out = write_csv(
    plastid_acc_data, 
    file_out("results/releases/ftol_plastid_accs.csv")),
  
  # - make table of gene partitions (start and end positions in alignment)
  plastid_gene_part_data = make_gene_part_table(plastid_genes_aligned_trimmed_renamed),
  
  # - write out table of gene partitions
  plastid_genes_part_data_out = write_csv(
    plastid_gene_part_data,
    file_out("results/releases/ftol_plastid_parts.csv")
  ),
  
  # - render data release README
  ftol_readme = rmarkdown::render(
    knitr_in("reports/results_readme/results_readme.Rmd"),
    output_file = here::here("results/releases/README.md"),
    quiet = TRUE
  ),

  # Generate reports ----

  plastid_tree_report = rmarkdown::render(
    knitr_in("reports/plastid_tree/plastid_tree.Rmd"),
    output_file = "plastid_tree.md",
    quiet = TRUE)

)
