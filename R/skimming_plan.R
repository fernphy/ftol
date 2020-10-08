# drake workflow plans

# Load data ----
data_plan <- drake_plan (
  
  # Read in names of Wei et al. 2017 protein-coding genes (83 genes).
  wei_genes_path = target("data_raw/wei_2017_coding_genes.txt", format = "file"),
  wei_genes = read_lines(wei_genes_path),
  
  # Read in PPGI taxonomic system.
  ppgi_taxonomy_path = target("data_raw/ppgi_taxonomy_mod.csv", format = "file"),
  ppgi_taxonomy = read_csv(ppgi_taxonomy_path),
  
  # Read in GoFlag metadata
  goflag_meta_path = target("data_raw/goflag/Pilot_Ferns_TargetCapture_Skimming.txt", format = "file"),
  goflag_meta = read_tsv(goflag_meta_path) %>%
    select(taxon = Taxon, targeted_capture_id = `Targeted Capture ID`, genome_skimming_id  = `Genome Skimming ID`) %>%
    filter(!is.na(taxon)),
  
  # extract vector of GoFlag sample IDs
  goflag_ids = goflag_meta %>%
    pivot_longer(names_to = "id_type", values_to = "id", -taxon) %>%
    pull(id)
  
)

# 01_trimmomatic ----

# Trim raw fastq files with trimmomatic.
trimmomatic_plan <- drake_plan (
  
  # Read in vector of raw fastq.gz files with absolute paths and
  # make table of arguments to feed into trimmomatic.
  raw_reads_dir = target("data_raw/goflag/", format = "file"),
  
  fastq_files = get_fastq_files_for_trimmomatic(
    fastq_raw_dir = raw_reads_dir,
    fastq_pattern = "fastq.gz"
  ),
  
  trimmomatic_arg_table = make_trimmomatic_arg_table(
    fastq_files = fastq_files,
    ids = goflag_ids,
    outdir = here::here("intermediates/trimmomatic/")
  ),
  
  # Run trimmomatic on all demultiplexed raw reads, 
  # save the raw stdout and stderr to a tibble.
  trimmomatic_results_raw = target(
    trimmomatic_pe(
      inputFile1 = trimmomatic_arg_table$inputFile1,
      inputFile2 = trimmomatic_arg_table$inputFile2,
      outputFile1P = trimmomatic_arg_table$outputFile1P,
      outputFile1U = trimmomatic_arg_table$outputFile1U,
      outputFile2P = trimmomatic_arg_table$outputFile2P,
      outputFile2U = trimmomatic_arg_table$outputFile2U,
      trim_settings = glue::glue(
        "ILLUMINACLIP:{here::here('data_raw/fasta/TruSeq3-PE-2.fa')}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50"
      ),
      threads = 1,
      adapters = file_in("data_raw/fasta/TruSeq3-PE-2.fa")
    ),
    dynamic = map(trimmomatic_arg_table
    )
  ),
  
  trimmomatic_results_summary = target(
    bind_rows(trimmomatic_results_raw),
    dynamic = group(trimmomatic_results_raw)
  )
  
)

# 02_format_targets ----

format_plastid_targets_plan <- drake_plan(
  
  # Assemble set of coding genes from GenBank plastome data
  # using map(accessions, ~fetch_genes_from_plastome(., wei_genes))
  # Run this command locally as juno seems to have a hard time keeping a connection
  plastid_targets_path = target("temp/plastid_targets.RDS", format = "file"),
  plastid_targets = readRDS(plastid_targets_path),
  
  # Collapse amino acid targets into single list, write it out
  plastid_aa_targets_out = ape::write.FASTA(
    transpose(plastid_targets)[["aa"]] %>% do.call(c, .),
    file_out("intermediates/hybpiper/plastid_aa_targets.fasta")
  ),
  
  # Write out plastid DNA targets as one fasta file per gene
  plastid_dna_targets_out = write_dna_targets_by_gene(
    plastid_targets, 
    here::here("intermediates/plastid_genes_ref")
  )
  
)

# 03_hybpiper ----
#
# Run hybpiper.
# Hybpiper will produce one folder for each sample in 03_hybpiper.

hybpiper_plan <- drake_plan (
  
  # Get vector of trimmed reads for hybpiper
  forward_reads = get_reads(
    data_dir = here::here("intermediates/trimmomatic/"),
    pattern = "R1.fastq",
    depends = trimmomatic_results_summary
  ),
  
  reverse_reads = get_reads(
    data_dir = here::here("intermediates/trimmomatic/"),
    pattern = "R2.fastq",
    depends = trimmomatic_results_summary
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
  plastid_hybpiper_results_each = target(
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
      bwa = FALSE),
    dynamic = map(paired_reads_list)
  ),
  
  # Combine results
  plastid_hybpiper_results = target(
    c(plastid_hybpiper_results_each),
    dynamic = group(plastid_hybpiper_results_each)
  ),
  
  plastid_samples = make_hybpiper_sample_file(
    in_dir = here::here("intermediates/hybpiper"), 
    pattern = "UFL|UFG", 
    out_path = file_out("intermediates/hybpiper/plastid_samples.txt"),
    depends = plastid_hybpiper_results
  ),
  
  plastid_lengths = get_seq_lengths(
    baitfile = here::here("intermediates/hybpiper/plastid_aa_targets.fasta"), 
    namelistfile = file_in("intermediates/hybpiper/plastid_samples.txt") %>% here::here(), 
    sequenceType = "aa",
    out_path = file_out("intermediates/hybpiper/plastid_lengths.txt"),
    wd = here::here("intermediates/hybpiper")
  ),
  
  plastid_stats = hybpiper_stats(
    seq_lengths = file_in("intermediates/hybpiper/plastid_lengths.txt") %>% here::here(), 
    namelistfile = file_in("intermediates/hybpiper/plastid_samples.txt") %>% here::here(),
    wd = here::here("intermediates/hybpiper")
  ),
  
  plastid_genes = retrieve_sequences(
    wd = here::here("intermediates/hybpiper/genes_recovered"),
    baitfile = here::here("intermediates/hybpiper/plastid_aa_targets.fasta"),
    sequence_dir = here::here("intermediates/hybpiper"), 
    sequenceType = "dna",
    depends = plastid_hybpiper_results),
  
  plastid_read_stats = get_read_stats(
    hybpiper_dir = here::here("intermediates/hybpiper"),
    depends = plastid_hybpiper_results
  ),
  
  # Align read fragments to reference
  plastid_read_fragments_aligned = align_hybpiper_reads_to_ref(
    hybpiper_dir = here::here("intermediates/hybpiper"),
    ref_dir = here::here("intermediates/plastid_genes_ref"),
    depends1 = plastid_hybpiper_results,
    depends2 = plastid_dna_targets_out
  ),
  
  # Extract consensus short reads (also excludes outliers)
  consensus_short_reads = extract_extract_short_reads_consensus_from_hybpiper(
    plastid_lengths = plastid_lengths,
    plastid_read_fragments_aligned = plastid_read_fragments_aligned)
  
)

# Combine plans ----
skimming_plan <- bind_plans(
  data_plan, 
  trimmomatic_plan,
  format_plastid_targets_plan,
  hybpiper_plan)
