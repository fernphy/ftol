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

tar_plan(
  # Prep for assembling Sanger plastid regions ----
  # Define variables used in plan:
  # - Target plastic loci (coding genes and spacers)
  target_loci = c(
    "atpA", "atpB", "matK", "rbcL", "rps4",
    "trnL-trnF", "rps4-trnS"),
  # - Most recent date cutoff for sampling genes
  date_cutoff = "2021/10/26",

  # Assemble initial reference Sanger sequences ----
  # These will be used for extracting target regions with BLAST.
  #
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
  # Trim sequences (lightly)
  fern_ref_seqs_trimmed = trim_genes(fern_ref_seqs_aligned),
  # Write out alignments for inspection FIXME
  tar_file(
    fern_ref_seqs_trimmed_out,
    write_fasta_tar(
      fern_ref_seqs_trimmed_dnabin,
      paste0("intermediates/ref_seqs/", target_loci, "_ref_aln.fasta")
    ),
    pattern = map(fern_ref_seqs_trimmed_dnabin, target_loci)
  ),

  # Download and extract Sanger sequences ----
  # Download raw fasta files
  tar_target(
    fern_sanger_seqs_raw,
    fetch_fern_sanger_seqs(
      target_loci, end_date = date_cutoff, accs_exclude_list = NULL),
    pattern = map(target_loci),
    # don't run in parallel, or will get HTTP status 429 errors
    deployment = "main"
  ),
  # Extract target regions
  tar_target(
    fern_sanger_extract_res,
    extract_from_ref_blast(
      query_seqtbl = fern_sanger_seqs_raw, # query: seqtbl with one row per sequence
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
    fetch_fern_metadata(target_loci, end_date = date_cutoff),
    pattern = map(target_loci),
    deployment = "main"
  ),
  raw_meta = unique(raw_meta_all),

  # Align sequences, second round ----
  # Filter to one best (longest) sequence per genus per target locus
  rep_fasta = filter_raw_fasta_by_genus(raw_fasta, raw_meta),
  # Align each target locus
  tar_target(
    fern_ref_seqs_aligned_2,
    align_ref_seqs(rep_fasta, target_loci, n_threads = 10),
    pattern = map(target_loci)
  ),
  # Trim sequences (lightly)
  # - dnabin format for IQTREE
  tar_target(
    fern_ref_seqs_trimmed_dnabin_2,
    trimal(
      fern_ref_seqs_aligned_2,
      c("-gt", "0.05", "-terminalonly"),
      return_seqtbl = FALSE),
    pattern = map(fern_ref_seqs_aligned_2)
  ),
  # - seqtbl format for downstream steps
  tar_target(
    fern_ref_seqs_trimmed_2,
    dnabin_to_seqtbl(fern_ref_seqs_trimmed_dnabin_2),
    pattern = map(fern_ref_seqs_trimmed_dnabin_2)
  ),
  # Build trees (to check quality of reference sequences)
  tar_target(
    fern_ref_seqs_tree_2,
    jntools::iqtree(
      fern_ref_seqs_trimmed_dnabin_2,
      m = "GTR+I+G", bb = 1000, nt = "AUTO",
      redo = TRUE, echo = TRUE, wd = tempdir()
      ),
    pattern = map(fern_ref_seqs_trimmed_dnabin_2)
  ),
  # Write out alignments for inspection
  tar_file(
    fern_ref_seqs_trimmed_out_2,
    write_fasta_tar(
      fern_ref_seqs_trimmed_dnabin_2,
      paste0("intermediates/ref_seqs/", target_loci, "_ref_aln_2.fasta")
    ),
    pattern = map(fern_ref_seqs_trimmed_dnabin_2, target_loci)
  ),
  # Write out trees for inspection
  tar_file(
    fern_ref_seqs_tree_out_2,
    write_tree_tar(
      fern_ref_seqs_tree_2,
      paste0("intermediates/ref_seqs/", target_loci, "_ref_phy_2.tree")
    ),
    pattern = map(fern_ref_seqs_tree_2, target_loci)
  )

)