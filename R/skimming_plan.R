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

)

# 01_fastp ----

# Trim raw fastq files with fastp
trimming_plan <- drake_plan (
  
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
      bwa = FALSE,
      # don't run exonerate or try to assemble genes
      other_args = c("--no-exonerate", "--no-assemble")),
    dynamic = map(paired_reads_list)
  ),
  
  # Combine results
  plastid_hybpiper_results = target(
    c(plastid_hybpiper_results_each),
    dynamic = group(plastid_hybpiper_results_each)
  ),
  
  # Align read fragments to reference
  # - loop over the list of samples
  each_extracted_reads_consensus = target(
    get_hybpip_consensus(
      sample = seq_cap_samples, 
      plastid_targets = plastid_targets,
      depends = plastid_hybpiper_results),
    dynamic = map(seq_cap_samples)
  ),
  
  # - combine results
  extracted_reads_consensus = target(
    bind_rows(each_extracted_reads_consensus),
    dynamic = group(each_extracted_reads_consensus)
  )
  
)

# Combine plans ----
skimming_plan <- bind_plans(
  data_plan, 
  trimming_plan,
  format_plastid_targets_plan,
  hybpiper_plan
  )
