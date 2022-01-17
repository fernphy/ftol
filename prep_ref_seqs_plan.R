library(targets)
library(tarchetypes)
source("R/packages.R")
source("R/functions.R")

# Specify path to folder with raw data
data_raw <- "_targets/user/data_raw"

# Set parallel back-end
plan(callr)

# Use targets workspaces for debugging
tar_option_set(workspace_on_error = TRUE)

tar_plan(
  # Sanger: Prep for assembling plastid regions ----
  # Define variables used in plan:
  # - Target plastic loci (coding genes and spacers)
  target_loci = c(
    "atpA", "atpB", "matK", "rbcL", "rps4",
    "trnL-trnF", "rps4-trnS"),
  # - Most recent date cutoff for sampling genes
  date_cutoff = "2021/12/31",
  # - Manually curated list of GenBank accessions to exclude from analysis
  tar_file(
    accs_exclude_path,
    path(data_raw, "accs_exclude.csv")),
  accs_exclude = read_csv(accs_exclude_path),

  # Sanger: Assemble initial reference sequences by parsing flatfiles ----
  # These will be used for extracting target regions with superCRUNCH

  # Download raw gene sequences for each target locus from GenBank
  tar_target(
    fern_ref_seqs_raw,
    fetch_fern_ref_seqs(
      target = target_loci,
      end_date = date_cutoff,
      # Only require intron for trnL-trnF
      req_intron = str_detect(target_loci, "trnL-trnF"),
      strict = FALSE,
      accs_exclude = accs_exclude$accession,
      workers = 48),
    pattern = map(target_loci),
    deployment = "main"
  ),
  # Filter to one longest seq per genus per target
  fern_ref_seqs = filter_ref_seqs(fern_ref_seqs_raw),
  # Align sequences
  tar_group_by(
    fern_ref_seqs_grouped,
    fern_ref_seqs,
    target # meaning target gene region, not a {targets} target
  ),
  tar_target(
    fern_ref_seqs_aligned,
    align_seqs_tbl(fern_ref_seqs_grouped),
    pattern = map(fern_ref_seqs_grouped)
  ),
  # Trim sequences
  fern_ref_seqs_trimmed_full_rps4 = trim_genes(fern_ref_seqs_aligned),
  # Trim out rps4 gene (only keep rps4-trnS spacer)
  fern_ref_seqs_trimmed = trim_align_by_motif(fern_ref_seqs_trimmed_full_rps4),
  # Write out alignments for inspection
  tar_file(
    fern_ref_seqs_trimmed_out,
    write_fasta_from_tbl(
      fern_ref_seqs_trimmed,
      dir = "intermediates/ref_seqs",
      prefix = "ref_aln_init_",
      postfix = ".fasta"),
    pattern = map(fern_ref_seqs_trimmed)
  ),

  # Sanger: Download and extract sequences with superCRUNCH ----

  # Download raw fasta files
  tar_target(
    fern_sanger_seqs_raw,
    fetch_fern_sanger_seqs(
      target_loci,
      end_date = date_cutoff,
      accs_exclude = accs_exclude$accession,
      strict = FALSE),
    pattern = map(target_loci),
    # don't run in parallel, or will get HTTP status 429 errors
    deployment = "main"
  ),
  # Extract target regions
  tar_target(
    fern_sanger_extract_res,
    extract_from_ref_blast(
      query_seqtbl = fern_sanger_seqs_raw, # query: seqtbl with one row per sequence # nolint
      ref_seqtbl = fern_ref_seqs_trimmed, # ref: tbl with one row per alignment
      target = target_loci,
      blast_flavor = "dc-megablast",
      other_args = c("-m", "span", "--threads", "4")
    ),
    pattern = map(target_loci)
  ),
  raw_fasta = clean_extract_res(fern_sanger_extract_res, "dc-megablast"),
  # Fetch metadata
  tar_target(
    raw_meta_all,
    fetch_fern_metadata(
      target_loci,
      end_date = date_cutoff,
      accs_exclude = accs_exclude$accession,
      strict = FALSE),
    pattern = map(target_loci),
    deployment = "main"
  ),
  raw_meta = unique(raw_meta_all),

  # Sanger: Clean final set of reference sequences ----

  # Filter to one best (longest) sequence per genus per target locus
  fern_ref_seqs_2 = filter_raw_fasta_by_genus(raw_fasta, raw_meta),
  # Align each target locus
  tar_group_by(
    fern_ref_seqs_grouped_2,
    fern_ref_seqs_2,
    target # meaning target gene region, not a {targets} target
  ),
  tar_target(
    fern_ref_seqs_aligned_2,
    align_seqs_tbl(fern_ref_seqs_grouped_2),
    pattern = map(fern_ref_seqs_grouped_2)
  ),
  # Rename sequences in alignment and trim genes
  fern_ref_seqs_trimmed_2 = fern_ref_seqs_aligned_2 %>%
    mutate(species = glue("{species}_{accession}")) %>%
    trim_genes(),
  # Write out final reference sequence alignments to raw data folder
  # to will use in main _targets.R plan
  tar_file(
    fern_ref_seqs_trimmed_out_2,
    write_fasta_from_tbl(
      fern_ref_seqs_trimmed_2,
      dir = fs::path(data_raw, "ref_aln"),
      prefix = "",
      postfix = ".fasta"),
    pattern = map(fern_ref_seqs_trimmed_2)
  ),
  # Build trees (to check quality of reference sequences)
  tar_target(
    fern_ref_seqs_tree_2,
    build_tree_from_alignment_df(
      fern_ref_seqs_trimmed_2,
      program = "iqtree"
    ),
    pattern = map(fern_ref_seqs_trimmed_2)
  ),
  # Write out trees for inspection
  tar_file(
    fern_ref_seqs_tree_out_2,
    write_tree_from_tbl(fern_ref_seqs_tree_2,
      dir = "intermediates/ref_seqs",
      prefix = "ref_phy_",
      postfix = ".tree"),
    pattern = map(fern_ref_seqs_tree_2)
  ),

  # Plastome: Prep for assembling plastome genes ----

  # Load outgroup plastome accessions
  tar_file(
    plastome_outgroups_path,
    path(data_raw, "plastome_outgroups.csv")),
  plastome_outgroups = read_csv(plastome_outgroups_path),
  # Load list of coding genes to extract from plastomes
  # (based on genes of Wei et al 2017, then trimmed to non-duplicated genes)
  tar_file(
    target_plastome_genes_path,
    path(data_raw, "target_coding_genes.txt")),
  target_plastome_genes_all = read_lines(target_plastome_genes_path),
  # Exclude target Sanger genes
  # (since already made above in Sanger workflow)
  target_plastome_genes = target_plastome_genes_all[
    !target_plastome_genes_all %in% target_loci
  ],
  # Load pteridocat taxonomic reference
  # FIXME: temporary work-around for loading pteridocat data
  # until {pteridocat} package is live
  tar_file(pteridocat_file, "working/pteridocat_2022-01-07.csv"),
  pteridocat = read_csv(pteridocat_file),
  # Parse reference names
  pc_ref_names = ts_parse_names(unique(pteridocat$scientificName)),

  # Plastome: extract genes to use as reference ----

  # Download metadata (accession numbers and species names) from GenBank
  plastome_metadata_raw = download_plastome_metadata(
    end_date = date_cutoff,
    outgroups = plastome_outgroups,
    strict = TRUE,
    accs_exclude = accs_exclude$accession),
  # Resolve species names in metadata
  # (will drop accession if name could not be resolved)
  plastome_metadata_renamed = resolve_pterido_plastome_names(
    plastome_metadata_raw, plastome_outgroups, pc_ref_names, pteridocat
  ),
  # Download plastome gene sequences from GenBank
  target_plastome_accessions = unique(plastome_metadata_renamed$accession),
  tar_target(
    plastome_genes_raw,
    fetch_fern_genes_from_plastome(
      genes = target_plastome_genes,
      accession = target_plastome_accessions),
    pattern = map(target_plastome_accessions),
    deployment = "main"),
  # Rename accessions (to 'species_accession'),
  # filter to one representative sequence per genus,
  # group by gene
  tar_group_by(
    plastome_genes_unaligned,
    rename_and_filter_raw_plastome_seqs(
      plastome_genes_raw,
      plastome_metadata_renamed
      ),
    target # here, 'target' refers to the gene
  ),
  # Align sequences by gene
  tar_target(
    plastome_genes_aligned,
    align_seqs_tbl(plastome_genes_unaligned),
    pattern = map(plastome_genes_unaligned)
  ),
  # Trim alignments
  plastome_genes_aligned_trimmed = trim_plastome_genes(
    plastome_genes_aligned),
  # Write out final reference sequence alignments to raw data folder
  # to will use in main _targets.R plan
  tar_file(
    plastome_genes_ref_out,
    write_fasta_from_tbl(
      plastome_genes_aligned_trimmed,
      dir = fs::path(data_raw, "ref_aln"),
      prefix = "",
      postfix = ".fasta"),
    pattern = map(plastome_genes_aligned_trimmed)
  )
)