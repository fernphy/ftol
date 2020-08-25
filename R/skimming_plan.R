# drake workflow plans

# Load data ----
data_plan <- drake_plan (
  
  # Read in names of Wei et al. 2017 protein-coding genes (83 genes).
  wei_genes_path = target("data_raw/wei_2017_coding_genes.txt", format = "file"),
  wei_genes = read_lines(wei_genes_path),
  
  # Read in Wei et al. 2017 plastome accessions.
  wei_accessions_path = target("data/wei_2017_accessions.txt", format = "file"),
  wei_accessions = read_lines(wei_accessions_path),
  
  # Read in PPGI taxonomic system.
  ppgi_taxonomy_path = target("data_raw/ppgi_taxonomy_mod.csv", format = "file"),
  ppgi_taxonomy = read_csv(ppgi_taxonomy_path),
  
  # Read in GoFlag metadata
  goflag_meta_path = target("data_raw/goflag/Pilot_Ferns_TargetCapture_Skimming.txt", format = "file"),
  goflag_meta = read_tsv(goflag_meta_path) %>%
    janitor::clean_names() %>%
    filter(!is.na(taxon)),
  
  # extract vector of GoFlag sample IDs
  goflag_ids = goflag_meta %>%
    select(taxon, targeted_capture_id, genome_skimming_id) %>%
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
        "ILLUMINACLIP:{here('data_raw/fasta/TruSeq3-PE-2.fa')}:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50"
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
  # Assemble set of 83 coding genes from GenBank plastome data
  # This part is memory intensive: requires ca. 60 gb per CPU
  wei_gene_list = gbfetch::assemble_gene_set(wei_accessions, wei_genes),
  
  # Make target file of 83 coding plastid genes
  plastid_dna_targets = make_plastid_target_file(
    gene_list = wei_gene_list,
    accessions = wei_accessions,
    gene_names = wei_genes,
    taxonomy_data = ppgi_taxonomy,
    out_path = file_out("intermediates/hybpiper/plastid_dna_targets.fasta")
  )
  
)

# 03_hybpiper ----
#
# Run hybpiper.
#
# Run in parallel with one CPU per sample, should take ca. 1 hr per sample.
# So if we run with 8 CPUs, will take ca. 6 hours.
# Hybpiper will produce one folder for each sample in 03_hybpiper.

hybpiper_plan <- drake_plan (
  
  # Get vector of trimmed reads for hybpiper
  forward_reads = get_reads(
    data_dir = here("intermediates/trimmomatic/"),
    pattern = "R1.fastq",
    depends = trimmomatic_results_summary
  ),
  
  reverse_reads = get_reads(
    data_dir = here("intermediates/trimmomatic/"),
    pattern = "R2.fastq",
    depends = trimmomatic_results_summary
  ),
  
  # Make list of paired reads for HybPiper
  paired_reads_list = make_paired_reads_list (
    forward_reads,
    reverse_reads
  ),
  
  # Plastid genes (85 coding genes adapted from Wei et al 2017)
  
  plastid_hybpiper_results_each = target(
    reads_first(
      wd = here("intermediates/hybpiper"),
      echo = FALSE,
      baitfile = file_in("intermediates/hybpiper/plastid_dna_targets.fasta"),
      # When paired_reads_list gets split up by dynamic mapping, it is split into lists.
      # The character vector we want for `readfiles` is the first element of each list
      readfiles = paired_reads_list[[1]],
      cpu = 1, 
      bwa = TRUE),
    dynamic = map(paired_reads_list)
  ),
  
  plastid_hybpiper_results = target(
    c(plastid_hybpiper_results_each),
    dynamic = group(plastid_hybpiper_results_each)
  ),
  
  plastid_samples = make_hybpiper_sample_file(
    in_dir = here("intermediates/hybpiper"), 
    pattern = "4938|JNG", 
    out_path = file_out("intermediates/hybpiper/plastid_samples.txt"),
    depends = plastid_hybpiper_results
  ),
  
  plastid_lengths = get_seq_lengths(
    baitfile = here("intermediates/hybpiper/plastid_dna_targets.fasta"), 
    namelistfile = file_in("intermediates/hybpiper/plastid_samples.txt"), 
    sequenceType = "dna",
    out_path = file_out("intermediates/hybpiper/plastid_lengths.txt"),
    wd = here("intermediates/hybpiper")
  ),
  
  plastid_stats = hybpiper_stats(
    seq_lengths = file_in("intermediates/hybpiper/plastid_lengths.txt"), 
    namelistfile = file_in("intermediates/hybpiper/plastid_samples.txt"),
    wd = here("intermediates/hybpiper")
  ),
  
  plastid_genes = retrieve_sequences(
    wd = here("intermediates/hybpiper/genes_recovered/"),
    baitfile = here("intermediates/hybpiper/plastid_dna_targets.fasta"),
    sequence_dir = here("intermediates/hybpiper"), 
    sequenceType = "dna",
    depends = plastid_hybpiper_results),
  
)

# Combine plans ----
skimming_plan <- bind_plans(
  data_plan, 
  trimmomatic_plan,
  format_plastid_targets_plan,
  hybpiper_plan)
