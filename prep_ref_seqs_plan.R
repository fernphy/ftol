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
      query_seqtbl = fern_sanger_seqs_raw,
      ref_seqtbl = fern_ref_seqs_trimmed_clean,
      target = target_loci,
      blast_flavor = "dc-megablast",
      other_args = c("-m", "span", "--threads", "4")
    ),
    pattern = map(target_loci)
  ),
  # after comparing results between blastn and
  # dc-megablast, dc-megablast seems to work better
  raw_fasta = clean_extract_res(fern_sanger_extract_res, "dc-megablast"),
  # Fetch metadata
  tar_target(
    raw_meta_all,
    fetch_fern_metadata(target_loci, end_date = date_cutoff),
    pattern = map(target_loci),
    deployment = "main"
  ),
  raw_meta = unique(raw_meta_all),

  # Align sequences ----
  rep_fasta = filter_raw_fasta_by_genus(raw_fasta, raw_meta)

)