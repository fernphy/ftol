# Plastome backbone investigation: Cyatheales + Salviniales topology

## Problem statement

Across many iterations of the plastome (whole-plastid-genome) backbone
analysis, the resulting tree kept recovering tree ferns (Cyatheales) as
sister to water ferns (Salviniales). This contradicts other published fern
phylogenies, which do not recover this relationship. Starting at commit
`b7e021e` ("Used partitioned analysis for plastome tree"), a series of
increasingly involved approaches were tried to see whether the placement was
an artifact of the analysis method (model choice, partitioning, alignment
quality, saturation at 3rd codon positions, insufficient/ misleading
backbone constraint) that could be corrected.

After exhausting the approaches below, none resolved the topology. The
working conclusion is that this may be a genuine signal in the plastome data
rather than an analysis artifact — a question that could become its own
focused paper later, taking a deep dive on this node specifically. This
document exists so that a future investigation doesn't have to reconstruct
what was tried from git history: every approach is listed chronologically
below with its rationale (where recorded) and the full code involved.

Following this investigation, the codebase was reverted to the
pre-`b7e021e`, non-partitioned analysis method (see commit reverting
`_targets.R`/`R/functions.R`) so the pipeline could proceed with a new
GenBank release. Nothing below reflects the current state of the code —
it's a historical record of abandoned attempts.

## Approaches tried (chronological)

### 1. Partitioned IQ-TREE analysis (`b7e021e`, 2026-06-22)

No commit body was recorded; rationale is inferred from the code alone.
Introduced one partition per locus (79 loci for the plastome dataset) so IQ-TREE
could fit a separate best-fit model per locus rather than one model for the
whole concatenated alignment, on the theory that per-locus model
misspecification could be contributing to the odd placement. Used
`m = "MFP+MERGE"` (fit models, then greedily merge partitions by BIC) and
constrained model search to `-mrate E,I,G,I+G` (no FreeRate).

```r
# R/functions.R — new function
#' Write an IQ-TREE partition file in RAxML format
#'
#' @param parts_table Tibble with columns "locus", "start", "end"
#' @param out_path Path to write the partition file
#'
#' @return Path to the written partition file
#'
write_iqtree_partition_file <- function(parts_table, out_path) {
  lines <- paste0(
    "DNA, ", parts_table$locus, " = ", parts_table$start, "-", parts_table$end
  )
  writeLines(lines, out_path)
  return(out_path)
}
```

```r
# R/functions.R — iqtree() wrapper: -spp flag renamed to -p
# (old: if (!is.null(spp)) "-spp")
if (!is.null(spp)) "-p",
fs::path_abs(spp),
```

```r
# _targets.R — partition file + partitioned plastome_tree
tar_file(
  plastome_partition_file,
  write_iqtree_partition_file(
    plastome_parts_table,
    path(int_dir, "iqtree/plastome/plastome_partitions.txt")
  )
),
tar_target(
  plastome_tree,
  iqtree(
    plastome_alignment,
    m = "MFP+MERGE", # merge partitions by BIC after model testing
    bb = 1000,
    nt = plastome_tree_nt_setting,
    seed = 20220123,
    redo = plastome_tree_redo_setting,
    echo = TRUE,
    wd = path(int_dir, "iqtree/plastome"),
    spp = plastome_partition_file,
    other_args = c(
      "-mset", "GTR", # only test GTR models
      "-mrate", "E,I,G,I+G", # don't test free-rate models
      "-t", "PARS",
      "--rcluster-max", "10" # limit partition merge search to top 10% of candidates
    ),
    tree_path = path(
      int_dir, "iqtree/plastome/plastome_alignment.phy.contree"
    )
  ),
  deployment = "main"
),
```

### 2. FreeRate models + apply partitioning to Sanger tree too (`ef6d48b`, 2026-06-22)

Full commit message (quoted directly — this is the most detailed rationale
recorded in the whole saga):

> Add one partition per locus for both the plastome (79 loci) and Sanger ML
> (7 loci) analyses so each locus gets its own best-fit substitution model.
> Remove the -mrate E,I,G,I+G constraint that was blocking FreeRate models;
> the previous constraint forced GTR+I+G4 on every partition, whereas recent
> literature finds GTR+R5 optimal for plastome data. Also fix the plastome
> tree_path to point at the partition-named IQ-TREE output
> (plastome_partitions.txt.contree) instead of the stale unpartitioned
> result.

```r
# _targets.R — plastome_tree: drop the -mrate constraint, fix tree_path
      spp = plastome_partition_file,
      other_args = c(
        "-mset", "GTR", # only test GTR models
        "-t", "PARS",
        "--rcluster-max", "10" # limit partition merge search to top 10% of candidates
      ),
      tree_path = path(
        int_dir, "iqtree/plastome/plastome_partitions.txt.contree"
      )
```

```r
# _targets.R — sanger_ml_tree_rep: drop -mrate constraint, add partition file
      m = "MFP", # run modelfinder and use best model
      other_args = c(
        "-mset", "GTR", # only test GTR family of models
        "-t", "PARS",
        "-g", path_abs(constraint_tree_file)
      ),
      bb = 1000,
      nt = sanger_ml_tree_nt_setting,
      seed = iqtree_sanger_seeds,
      redo = sanger_ml_tree_redo_setting,
      wd = iqtree_sanger_dirs,
      spp = sanger_partition_file,
      tree_path = c(
        ml_tree = path(iqtree_sanger_dirs, "sanger_partitions.txt.treefile"),
        con_tree = path(iqtree_sanger_dirs, "sanger_partitions.txt.contree")
      ),
```

```r
# _targets.R — new sanger partition file target
tar_file(
  sanger_partition_file,
  write_iqtree_partition_file(
    sanger_parts_table,
    path(int_dir, "iqtree/sanger_partitions.txt")
  )
),
```

### 3. IQ-TREE wrapper bug fix: `--redo` needs two dashes (`68f2363`, 2026-06-25)

Pure bug fix, not a topology attempt. Commit message: "'redo' needed two
dashes."

```r
# R/functions.R — iqtree()
    if (isTRUE(redo)) "--redo",   # was: if (isTRUE(redo)) "-redo",
```

### 4. Disable partition merging (`7108b04`, 2026-06-25)

Commit message: "Too computationally expensive" (this is the entire recorded
rationale). Reverted the `MFP+MERGE` model setting back to plain `MFP` and
removed `--rcluster-max 10`, so partitions from step 1 are still used but no
longer merged by BIC.

```r
# _targets.R — plastome_tree
      m = "MFP", # test model followed by ML analysis   (was "MFP+MERGE")
      ...
      other_args = c(
        "-mset", "GTR", # only test GTR models
        "-t", "PARS"
        # ("--rcluster-max", "10" removed)
      ),
```

### 5. First attempt at excluding 3rd codon position: `trim_cds_no3rd` (`cfa395f`, 2026-07-03)

No commit body. Rationale inferred from code and docstrings: 3rd codon
positions in protein-coding plastome loci are expected to be saturated with
substitutions at this phylogenetic depth, which can produce long-branch
attraction artifacts (a classic cause of exactly the kind of spurious
sister-group relationship being chased here). This first version simply
stripped every 3rd column from the existing (non-codon-aware) MAFFT
alignment, on the assumption the reading frame was already intact, with a
guard (`cds_in_frame()`) that skips 3rd-position removal for a locus if doing
so introduces internal stop codons.

```r
# R/functions.R — new functions
#' Assert that a CDS alignment is in reading frame 0
#'
#' Selects the reference sequence with the fewest gaps, removes gaps, trims to
#' a multiple of 3, translates, and asserts no internal stop codons.
#'
#' @param aln_mat DNAbin matrix (rows = sequences, cols = alignment positions)
#' @param locus Character; locus name used in error messages
#'
cds_in_frame <- function(aln_mat, locus) {
  aln_char <- as.character(aln_mat)
  n_gaps <- rowSums(aln_char == "-")
  ref_char <- aln_char[which.min(n_gaps), ]
  ref_ungapped <- ref_char[ref_char != "-"]
  nc <- length(ref_ungapped)
  if (nc < 3) return(TRUE)
  nc_trim <- nc - nc %% 3
  ref_mat <- matrix(ref_ungapped[seq_len(nc_trim)], nrow = 1)
  ref_dnabin <- ape::as.DNAbin(ref_mat)
  aa <- as.character(ape::trans(ref_dnabin))[1L, ]
  internal_stops <- sum(aa[-length(aa)] == "*")
  if (internal_stops > 0L) {
    warning(glue::glue(
      "{locus}: reference sequence has {internal_stops} internal stop codon(s) ",
      "in reading frame 0 — skipping 3rd-position removal for this locus"
    ))
    return(FALSE)
  }
  TRUE
}

#' Trim CDS alignments after removing third codon positions
#'
#' Removes every third column from each gene's MAFFT alignment (while the
#' reading frame is still intact), then trims with trimal. Intended for
#' CDS-only input (no spacers).
#'
#' @param plastid_aligned Tibble output of align_seqs_tbl(), with columns
#'   seq, species, target, accession
#' @param name_col_in Name of column to use as sequence labels
#'
#' @return Tibble with columns "target" and "align_trimmed" (DNAbin matrices)
#'
trim_cds_no3rd <- function(plastid_aligned, name_col_in = "species") {
  plastid_aligned %>%
    select(seq, species, target, accession) %>%
    group_by(target) %>%
    nest(data = c(seq, species, accession)) %>%
    mutate(
      align_trimmed = map2(
        data, target,
        function(d, locus) {
          aln <- seqtbl_to_dnabin(d, name_col = name_col_in, seq_col = "seq")
          aln <- as.matrix(aln)
          n <- ncol(aln)
          if (n %% 3 != 0) {
            warning(glue::glue(
              "{locus}: MAFFT alignment has {n} columns (not divisible by 3); ",
              "trimming {n %% 3} column(s) from end to restore reading frame"
            ))
            n <- n - (n %% 3)
            aln <- aln[, seq_len(n)]
          }
          if (cds_in_frame(aln, locus)) {
            aln <- aln[, setdiff(seq_len(n), seq(3L, n, by = 3L))]
          }
          trimal(aln, other_args = c("-gt", "0.05"), return_seqtbl = FALSE)
        }
      )
    ) %>%
    ungroup() %>%
    select(-data)
}
```

```r
# _targets.R — new upstream targets feeding into plastome_partition_file_no3rd
plastid_genes_trimmed_no3rd = trim_cds_no3rd(plastid_genes_aligned),
plastome_alignment_no3rd = concatenate_to_ape(plastid_genes_trimmed_no3rd),
plastome_parts_table_no3rd = make_parts_table(
  plastid_genes_trimmed_no3rd, plastome_alignment_no3rd
),
tar_file(
  plastome_partition_file_no3rd,
  write_iqtree_partition_file(
    plastome_parts_table_no3rd,
    path(int_dir, "iqtree/plastome/plastome_partitions_no3rd.txt")
  )
),
```

Note: this `plastome_partition_file_no3rd` target was created here but never
actually wired into a tree-building `iqtree()` call in this commit — the
no-3rd-position alignment was prepared but not yet used for a tree.

### 6. Codon-aware alignment, manual implementation: `codon_align_and_trim_no3rd` (`53f1b06`, 2026-07-03)

No commit body. Reason for superseding approach 5 isn't recorded, but the
implementation itself is telling: instead of assuming the existing per-locus
MAFFT nucleotide alignment was already in frame, this version regenerates
the alignment in amino-acid space (translate → align proteins with MAFFT →
back-translate to nucleotides) so the alignment itself respects codon
boundaries, rather than just removing columns from a nucleotide alignment
that might have frameshift-inducing indels.

```r
# R/functions.R — replaces trim_cds_no3rd/cds_in_frame as the no3rd pipeline
codon_align_and_trim_no3rd <- function(plastid_aligned, name_col_in = "species") {
  plastid_aligned %>%
    select(seq, species, target, accession) %>%
    group_by(target) %>%
    nest(data = c(seq, species, accession)) %>%
    mutate(
      align_trimmed = map2(
        data, target,
        function(d, locus) {
          # 1. Strip alignment gaps to recover original unaligned sequences.
          #    Use accession as identifier (no spaces, safe for FASTA headers).
          aln <- seqtbl_to_dnabin(d, name_col = "accession", seq_col = "seq")
          seqs_ungapped <- as.list(ape::del.gaps(as.matrix(aln)))
          accs <- names(seqs_ungapped)
          acc_to_sp <- setNames(d[[name_col_in]], d$accession)

          # 2. Extract nt characters and trim to codon boundary
          nt_chars <- lapply(accs, function(acc) {
            ch <- as.character(seqs_ungapped[[acc]])
            r <- length(ch) %% 3L
            if (r != 0L) {
              warning(glue::glue(
                "{locus}/{acc}: trimming {r} trailing nt(s) to codon boundary"
              ))
              ch <- ch[seq_len(length(ch) - r)]
            }
            ch
          })
          names(nt_chars) <- accs

          # 3. Translate to AA strings; drop trailing stop codon (normal CDS end)
          aa_strs <- lapply(accs, function(acc) {
            ch <- nt_chars[[acc]]
            if (length(ch) < 3L) return("")
            aa <- as.character(
              ape::trans(ape::as.DNAbin(matrix(ch, nrow = 1L)))
            )[1L, ]
            if (length(aa) > 0L && aa[length(aa)] == "*") aa <- aa[-length(aa)]
            paste(aa, collapse = "")
          })
          names(aa_strs) <- accs

          # 4. Write AA FASTA and align with MAFFT in protein mode
          aa_in <- tempfile(fileext = ".faa")
          writeLines(
            unlist(lapply(accs, function(acc) c(paste0(">", acc), aa_strs[[acc]]))),
            aa_in
          )
          aa_out <- tempfile(fileext = "_aln.faa")
          system2("/usr/bin/mafft", c("--amino", "--quiet", aa_in), stdout = aa_out)

          # 5. Parse aligned AA FASTA
          raw <- readLines(aa_out)
          hdr_idx <- which(startsWith(raw, ">"))
          aa_aln <- setNames(
            lapply(seq_along(hdr_idx), function(i) {
              s <- hdr_idx[i] + 1L
              e <- if (i < length(hdr_idx)) hdr_idx[i + 1L] - 1L else length(raw)
              strsplit(paste(raw[s:e], collapse = ""), "")[[1L]]
            }),
            sub("^>", "", raw[hdr_idx])
          )
          aln_len_aa <- length(aa_aln[[accs[1L]]])

          # 6. Back-translate: each gap AA → "---", j-th non-gap AA → 3 nt
          nt_mat <- matrix(
            "-", nrow = length(accs), ncol = aln_len_aa * 3L,
            dimnames = list(accs, NULL)
          )
          for (acc in accs) {
            aa_vec <- aa_aln[[acc]]
            nt_src <- nt_chars[[acc]]
            non_gap <- which(aa_vec != "-")
            for (j in seq_along(non_gap)) {
              col <- (non_gap[j] - 1L) * 3L + 1L
              src_s <- (j - 1L) * 3L + 1L
              src_e <- src_s + 2L
              if (src_e <= length(nt_src)) {
                nt_mat[acc, col:(col + 2L)] <- nt_src[src_s:src_e]
              }
            }
          }

          # 7. Relabel rows from accession to species
          rownames(nt_mat) <- acc_to_sp[rownames(nt_mat)]

          # 8. Remove 3rd codon positions (columns 3, 6, 9, ...)
          n_cols <- ncol(nt_mat)
          keep_cols <- setdiff(seq_len(n_cols), seq(3L, n_cols, by = 3L))
          nt_mat <- nt_mat[, keep_cols, drop = FALSE]

          # 9. Trim with trimal
          trimal(
            ape::as.DNAbin(nt_mat),
            other_args = c("-gt", "0.05"),
            return_seqtbl = FALSE
          )
        }
      )
    ) %>%
    ungroup() %>%
    select(-data)
}
```

```r
# _targets.R
plastid_genes_trimmed_no3rd = codon_align_and_trim_no3rd(plastid_genes_aligned),
```

### 7. Switch to DECIPHER for codon-aware alignment (`8171fc6`, 2026-07-03)

No commit body beyond the subject line. Replaced the hand-rolled
translate/MAFFT-align/back-translate pipeline from approach 6 with
`DECIPHER::AlignTranslation()`, which does the same conceptual operation
(codon-aware alignment) using a purpose-built Bioconductor tool instead of
manual FASTA round-tripping through an external MAFFT process — likely a
robustness/maintainability motivation rather than a response to a specific
topology result (not recorded either way). Split the pipeline into two
functions: `codon_align_seqs_tbl()` (alignment only) and `strip_3rd_and_trim()`
(3rd-position removal + trimal), so the pre-removal codon alignment could
also be inspected directly. This commit also bumped `DECIPHER`/`Biostrings`
and several other Bioconductor/CRAN package versions in `renv.lock` as a side
effect of adding the new dependency.

```r
# R/functions.R
codon_align_seqs_tbl <- function(plastid_aligned, name_col_in = "species") {
  plastid_aligned %>%
    select(seq, species, target, accession) %>%
    group_by(target) %>%
    nest(data = c(seq, species, accession)) %>%
    mutate(
      align_trimmed = map(
        data,
        function(d) {
          # Strip per-sequence gaps to recover original unaligned sequences
          aln_char <- as.character(
            as.matrix(seqtbl_to_dnabin(d, name_col = "accession", seq_col = "seq"))
          )
          seqs_str <- apply(aln_char, 1, function(row) {
            paste(row[row != "-"], collapse = "")
          })
          acc_to_sp <- setNames(d[[name_col_in]], d$accession)

          # Codon-aware alignment via DECIPHER (input must be uppercase).
          # readingFrame = NA (default): auto-detect per sequence.
          aln_ss <- DECIPHER::AlignTranslation(
            Biostrings::DNAStringSet(toupper(seqs_str))
          )

          # Convert aligned DNAStringSet → lowercase character matrix → DNAbin
          aln_mat <- do.call(rbind, strsplit(tolower(as.character(aln_ss)), ""))
          rownames(aln_mat) <- acc_to_sp[names(aln_ss)]
          ape::as.DNAbin(aln_mat)
        }
      )
    ) %>%
    ungroup() %>%
    select(-data)
}

strip_3rd_and_trim <- function(codon_aligned_tbl) {
  codon_aligned_tbl %>%
    mutate(
      align_trimmed = map(
        align_trimmed,
        function(aln) {
          n <- ncol(aln)
          keep <- setdiff(seq_len(n), seq(3L, n, by = 3L))
          trimal(aln[, keep, drop = FALSE],
                 other_args = c("-gt", "0.05"),
                 return_seqtbl = FALSE)
        }
      )
    )
}
```

```r
# _targets.R
# Step 1: codon-aware alignment per locus (full 3 positions, no trimal yet)
plastid_genes_aligned_codon = codon_align_seqs_tbl(plastid_genes_aligned),
# Step 2: concatenate for visual inspection
plastome_alignment_codon = concatenate_to_ape(plastid_genes_aligned_codon),
# Step 3: strip 3rd codon positions and trim
plastid_genes_trimmed_no3rd = strip_3rd_and_trim(plastid_genes_aligned_codon),
```

### 8. MAFFT `--adjustdirection` pass before DECIPHER (`17cdc3a`, 2026-07-03)

Commit message: "Use mafft before DECIPHER to make sure 5-3 direction is
correct." Refactored `codon_align_seqs_tbl()` into `codon_align_locus()`
(one call per locus, run as a branched `targets` pattern instead of an
internal loop) and, per its own docstring, added a preliminary
`ips::mafft(..., options = "--adjustdirection")` pass: sequences are aligned
per-family upstream, so orientation is consistent *within* a family but may
differ *across* families — this step brings all sequences for a locus to a
consensus strand direction before codon-aware alignment, since
`DECIPHER::AlignTranslation()` does not itself correct strand orientation.
Also added an unused helper, `sample_loci_fasta()`, for spot-checking a
random subset of loci from the codon alignment.

```r
# R/functions.R
#' Codon-aware alignment for a single plastid locus
#'
#' Strips alignment gaps from a per-family MAFFT alignment to recover unaligned
#' sequences, normalizes strand orientation across all families with a
#' preliminary MAFFT --adjustdirection pass, then performs codon-aware
#' alignment via DECIPHER::AlignTranslation(). Intended to be called as a
#' branched targets target, one branch per locus.
#'
#' @param locus_tbl Tibble for a single locus with columns seq, species,
#'   target, accession (one row per sequence)
#' @param name_col_in Name of column to use as sequence labels in output
#'
#' @return One-row tibble with columns "target" and "align_trimmed"
#'   (DNAbin matrix)
#'
codon_align_locus <- function(locus_tbl, name_col_in = "species") {
  locus <- unique(locus_tbl$target)
  acc_to_sp <- setNames(locus_tbl[[name_col_in]], locus_tbl$accession)

  # Strip per-sequence gaps to recover original unaligned sequences
  aln_char <- as.character(
    as.matrix(
      seqtbl_to_dnabin(locus_tbl, name_col = "accession", seq_col = "seq")
    )
  )
  seqs_str <- apply(aln_char, 1, function(row) {
    paste(row[row != "-"], collapse = "")
  })

  # Normalize orientation with MAFFT --adjustdirection.
  # plastid_genes_aligned is built per-family, so orientation is consistent
  # within each family group but may differ across families. This step
  # brings all sequences to a consensus direction before codon alignment.
  tmp_in <- tempfile(fileext = ".fasta")
  writeLines(
    unlist(lapply(names(seqs_str), function(nm) c(paste0(">", nm), seqs_str[[nm]]))),
    tmp_in
  )
  seqs_dnabin <- ape::read.FASTA(tmp_in)
  dir_aln <- ips::mafft(seqs_dnabin, options = "--adjustdirection",
                        exec = "/usr/bin/mafft")

  # Strip alignment gaps; drop _R_ suffix MAFFT adds to reversed sequences
  dir_char <- as.character(as.matrix(dir_aln))
  rownames(dir_char) <- str_remove_all(rownames(dir_char), "_R_")
  seqs_corrected <- apply(dir_char, 1, function(row) {
    paste(row[row != "-"], collapse = "")
  })

  # Codon-aware alignment via DECIPHER (input must be uppercase).
  # readingFrame = NA (default): auto-detect per sequence.
  aln_ss <- DECIPHER::AlignTranslation(
    Biostrings::DNAStringSet(toupper(seqs_corrected))
  )

  # Convert aligned DNAStringSet → lowercase character matrix → DNAbin
  aln_mat <- do.call(rbind, strsplit(tolower(as.character(aln_ss)), ""))
  rownames(aln_mat) <- acc_to_sp[names(aln_ss)]

  tibble::tibble(
    target = locus,
    align_trimmed = list(ape::as.DNAbin(aln_mat))
  )
}

#' Strip 3rd codon positions and trim a codon-aligned tibble
#'
#' Removes every third column (positions 3, 6, 9, ...) from each locus
#' alignment, then trims with trimal (gap threshold 0.05). Intended to be
#' applied to the output of codon_align_locus(), where codon boundaries are
#' guaranteed to be intact.
#'
#' @param codon_aligned_tbl Tibble with columns "target" and "align_trimmed"
#'   (DNAbin matrices from codon_align_locus())
#'
#' @return Tibble with the same structure, with 3rd positions removed and
#'   columns trimmed by trimal
#'
strip_3rd_and_trim <- function(codon_aligned_tbl) {
  # (body unchanged from approach 7)
}

#' Write a random subset of loci from a per-locus alignment tibble to FASTA
#'
#' Concatenates a random sample of loci and writes the result to a FASTA file.
#' Useful for spot-checking a codon-aware alignment without loading the full
#' concatenated matrix.
#'
#' @param aln_tbl Tibble with columns "target" and "align_trimmed" (one row
#'   per locus), e.g. the output of codon_align_locus()
#' @param out_path Path to write the FASTA file
#' @param n Number of loci to sample (default 20)
#' @param seed Optional random seed for reproducibility
#'
#' @return Path to the written FASTA file (invisibly)
#'
sample_loci_fasta <- function(aln_tbl, out_path, n = 20, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sampled <- aln_tbl[sample(nrow(aln_tbl), min(n, nrow(aln_tbl))), ]
  aln <- concatenate_to_ape(sampled)
  write_fasta_tar(aln, out_path)
}
```

```r
# _targets.R
tar_group_by(
  plastid_genes_aligned_by_locus,
  plastid_genes_aligned,
  target
),
tar_target(
  plastid_genes_aligned_codon,
  codon_align_locus(plastid_genes_aligned_by_locus),
  pattern = map(plastid_genes_aligned_by_locus)
),
plastome_alignment_codon = concatenate_to_ape(plastid_genes_aligned_codon),
tar_file(
  plastome_alignment_codon_file,
  write_fasta_tar(
    plastome_alignment_codon,
    path(int_dir, "plastome_alignment_codon.fasta")
  )
),
plastid_genes_trimmed_no3rd = strip_3rd_and_trim(plastid_genes_aligned_codon),
```

### 9. One-taxon-per-genus exemplar backbone constraint tree (`47aa196`, 2026-07-04)

No commit body. Rationale from the function docstring: rather than build the
backbone tree from every species (where within-genus noise or
misidentification could be contributing to the wrong placement), build a
reduced backbone tree from one exemplar species per genus (the
most-data-complete species in each multi-species genus, using the
no-3rd-position alignment), infer a tree from that smaller, presumably
cleaner dataset, then use *that* backbone tree as a topology constraint
(`-g`) for the full all-species plastome tree. Also added general
`g`/constraint-tree support and `wd` auto-creation to the `iqtree()` wrapper.

```r
# R/functions.R
#' Select exemplar species for the plastome backbone tree
#'
#' For each genus recognised by PPG with more than one species in the plastome
#' dataset, selects the single most data-complete species (most loci with any
#' sequence data) as an exemplar, unless the genus is listed in
#' \code{non_mono_genera}. All species from singleton genera, unrecognised
#' genera, and flagged non-monophyletic genera are retained in full.
#'
#' @param plastid_genes_trimmed Tibble with columns "target" and
#'   "align_trimmed" (DNAbin matrices), e.g. output of strip_3rd_and_trim()
#' @param ppgi_taxonomy Tibble of PPG genus-level taxonomy as produced by
#'   taxlist_to_df(), with a column "genus"
#' @param non_mono_genera Character vector of genus names known or suspected
#'   to be non-monophyletic; all species from these genera are retained
#'
#' @return Character vector of selected species names (sorted)
#'
select_plastome_exemplars <- function(
  plastid_genes_trimmed,
  ppgi_taxonomy,
  non_mono_genera = character(0)
) {
  # Count loci with any sequence data per species
  n_loci <- plastid_genes_trimmed$align_trimmed %>%
    lapply(function(mat) {
      aln_char <- as.character(mat)
      rownames(mat)[rowSums(aln_char != "-") > 0]
    }) %>%
    unlist() %>%
    table()

  sp_data <- tibble::tibble(
    species = names(n_loci),
    n_loci  = as.integer(n_loci),
    genus   = stringr::word(species, 1, sep = "_")
  )

  ppg_genera <- unique(ppgi_taxonomy$genus)

  # Genera with >1 species, recognised by PPG, not flagged as non-monophyletic
  mono_multi_genera <- sp_data %>%
    dplyr::group_by(genus) %>%
    dplyr::filter(
      dplyr::n() > 1,
      genus %in% ppg_genera,
      !genus %in% non_mono_genera
    ) %>%
    dplyr::pull(genus) %>%
    unique()

  exemplars <- sp_data %>%
    dplyr::filter(genus %in% mono_multi_genera) %>%
    dplyr::group_by(genus) %>%
    dplyr::slice_max(n_loci, n = 1, with_ties = FALSE) %>%
    dplyr::pull(species)

  others <- sp_data %>%
    dplyr::filter(!genus %in% mono_multi_genera) %>%
    dplyr::pull(species)

  sort(unique(c(exemplars, others)))
}
```

```r
# R/functions.R — iqtree() wrapper additions
iqtree <- function(
  ...,
  spp = NULL,
  g = NULL,               # new: path to constraint tree for -g
  ...
) {
  ...
  if (!is.null(g)) {
    assertthat::assert_that(assertthat::is.readable(g))
  }
  ...
  wd <- fs::path_norm(wd)
  fs::dir_create(wd, recurse = TRUE)   # new: auto-create working dir
  ...
  args <- c(
    ...,
    if (!is.null(g)) "-g",
    if (!is.null(g)) fs::path_abs(g),
    ...
  )
}
```

```r
# _targets.R
non_mono_genera_plastome <- character(0)

plastome_exemplars = select_plastome_exemplars(
  plastid_genes_trimmed_no3rd,
  ppgi_taxonomy,
  non_mono_genera = non_mono_genera_plastome
),
plastome_alignment_backbone = plastome_alignment_no3rd[
  rownames(plastome_alignment_no3rd) %in% plastome_exemplars, ,
  drop = FALSE
],
tar_file(
  plastome_partition_file_backbone,
  write_iqtree_partition_file(
    plastome_parts_table_no3rd,
    path(int_dir, "iqtree/plastome_backbone/plastome_partitions_no3rd.txt")
  ),
  deployment = "main"
),
tar_target(
  plastome_tree_backbone,
  iqtree(
    plastome_alignment_backbone,
    spp = plastome_partition_file_backbone,
    m = "MFP",
    bb = 1000,
    nt = plastome_tree_nt_setting,
    seed = 20220123,
    redo = plastome_backbone_redo_setting,
    echo = TRUE,
    wd = path(int_dir, "iqtree/plastome_backbone"),
    other_args = c("-t", "PARS"),
    tree_path = path(
      int_dir,
      "iqtree/plastome_backbone/plastome_partitions_no3rd.txt.contree"
    )
  ),
  deployment = "main"
),
tar_file(
  plastome_backbone_constraint_file,
  write_tree_tar(
    plastome_tree_backbone,
    path(int_dir, "iqtree/plastome_backbone/backbone.tre")
  )
),
# Full plastome tree: all species, with 3rd positions, constrained by backbone
tar_target(
  plastome_tree,
  iqtree(
    plastome_alignment,
    spp = plastome_partition_file,
    m = "MFP",
    ...
    g = plastome_backbone_constraint_file,
    other_args = c("-t", "PARS"),
    tree_path = path(
      int_dir, "iqtree/plastome/plastome_partitions.txt.contree"
    )
  ),
  deployment = "main"
),
```

### 10. Du et al. (2022) accession-matched backbone constraint (`05785d4`, 2026-07-04)

No commit body. Replaced the exemplar-selection approach (9) with a
different backbone-construction strategy: rather than picking one
representative per genus from *our own* dataset, the accessions used in Du
et al. (2022) — an independent, presumably well-vetted plastome phylogeny —
were mapped onto FTOL's GenBank-equivalent accessions (via
`data_raw/du2022_accessions.csv`, a manually curated NC_-to-GenBank mapping),
and the FTOL plastome alignment was subset to just those species. A backbone
tree was built from that subset using `GTR+F+R5` directly (skipping
ModelFinder) and used as the `-g` constraint for the full plastome tree, on
the theory that constraining to a topology independently derived by domain
experts elsewhere might correct any artifact specific to FTOL's own
backbone-tree construction.

```r
# R/functions.R — new function
filter_genes_to_species <- function(plastid_genes_aligned_codon, species) {
  plastid_genes_aligned_codon %>%
    dplyr::mutate(
      align_trimmed = purrr::map(
        align_trimmed,
        function(aln) {
          keep <- rownames(aln) %in% species
          aln[keep, , drop = FALSE]
        }
      )
    )
}
```

```r
# _targets.R
tar_file_read(
  du2022_accession_list,
  path(data_raw, "du2022_accessions.csv"),
  read_csv(!!.x, col_types = cols(.default = "c"))
),
...
# Du et al. (2022) reference tree used as backbone constraint
# Species: subset of FTOL plastome dataset matching Du et al. accessions
du2022_species = plastome_metadata_renamed %>%
  dplyr::filter(accession %in% du2022_accession_list$accession) %>%
  dplyr::pull(species),
# Filter codon-aligned loci to Du et al. species
du2022_genes = filter_genes_to_species(
  plastid_genes_aligned_codon,
  du2022_species
),
du2022_alignment = concatenate_to_ape(du2022_genes),
tar_target(
  du2022_tree,
  iqtree(
    du2022_alignment,
    m = "GTR+F+R5",
    bb = 1000,
    nt = plastome_tree_nt_setting,
    seed = 20220123,
    redo = du2022_redo_setting,
    echo = TRUE,
    wd = path(int_dir, "iqtree/du2022"),
    other_args = c("-t", "PARS"),
    tree_path = path(int_dir, "iqtree/du2022/du2022_alignment.phy.contree")
  ),
  deployment = "main"
),
tar_file(
  plastome_backbone_constraint_file,
  write_tree_tar(
    du2022_tree,
    path(int_dir, "iqtree/du2022/backbone.tre")
  )
),
```

Note: `select_plastome_exemplars()` (approach 9) was left in `R/functions.R`
unused after this commit — it was not called again.

### 11. Order-level constraint tree: Cyatheales+Polypodiales vs. Salviniales (uncommitted, abandoned)

This final attempt was never committed — it existed only as a working-tree
change at the time this document was written, and is preserved here since it
represents a real negative result. Rather than build any intermediate
backbone tree, this approach directly constrains the full plastome tree's
topology via a hand-written 3-part Newick constraint: every Cyatheales
species plus every Salviniales species (all of them, not just one
representative), split into two clades, with a single Polypodiales
representative forcing Cyatheales to group with eupolypods to the exclusion
of Salviniales. The function's own docstring records why this had to include
*every* species from both orders rather than one representative each: IQ-TREE's
`-g` constraint only constrains the taxa actually named in the constraint
tree — with only one representative per order, all the *other* Cyatheales
and Salviniales species remain unconstrained and are free to recover
whatever topology (including the artifact-suspected sister relationship)
the data/model would otherwise produce.

```r
# R/functions.R
#' Write an order-level constraint tree targeting Cyatheales placement
#'
#' Generates a Newick constraint for IQ-TREE (\code{-g}) that enforces the
#' known placement of Cyatheales and Salviniales: all Cyatheales form a clade
#' with eupolypods (Polypodiales) to the exclusion of all Salviniales, and all
#' Salviniales form their own monophyletic clade. Using one representative per
#' order leaves the remaining taxa free and does not work because IQ-TREE only
#' constrains the named species — all other Cyatheales/Salviniales remain free
#' to follow the long-branch artefact. Including every species from both orders
#' prevents that.
#'
#' @param alignment DNAbin matrix whose row names are species names (underscore-
#'   separated, e.g. \code{Alsophila_costularis}).
#' @param ppgi_taxonomy Data frame with columns \code{genus} and \code{order}
#'   (one row per genus, from PPG I taxonomy).
#' @param out_path File path for the output Newick file.
#'
#' @return \code{out_path} (invisibly).
#'
write_order_constraint_tree <- function(alignment, ppgi_taxonomy, out_path) {
  sp_names <- rownames(alignment)
  genera    <- stringr::word(sp_names, 1, sep = "_")

  sp_orders <- tibble::tibble(species = sp_names, genus = genera) %>%
    dplyr::left_join(
      dplyr::select(ppgi_taxonomy, genus, order),
      by = "genus"
    )

  get_all <- function(target_order) {
    sp_orders %>%
      dplyr::filter(order == target_order) %>%
      dplyr::pull(species)
  }

  pick_rep <- function(target_order) {
    sp_orders %>%
      dplyr::filter(order == target_order) %>%
      dplyr::slice(1) %>%
      dplyr::pull(species)
  }

  cyatheales_spp  <- get_all("Cyatheales")
  salviniales_spp <- get_all("Salviniales")
  polypodiales_sp <- pick_rep("Polypodiales")

  assertthat::assert_that(
    length(cyatheales_spp) > 0,
    msg = "No Cyatheales species found in alignment"
  )
  assertthat::assert_that(
    length(salviniales_spp) > 0,
    msg = "No Salviniales species found in alignment"
  )
  assertthat::assert_that(
    length(polypodiales_sp) == 1,
    msg = "No Polypodiales representative found in alignment"
  )

  cyath_str <- paste(cyatheales_spp,  collapse = ",")
  salv_str  <- paste(salviniales_spp, collapse = ",")
  newick    <- sprintf("((%s,%s),(%s));", cyath_str, polypodiales_sp, salv_str)

  fs::dir_create(fs::path_dir(out_path), recurse = TRUE)
  writeLines(newick, out_path)
  invisible(out_path)
}
```

```r
# _targets.R
# Order-level constraint tree targeting Cyatheales placement
# Minimal 3-leaf Newick: (Cyatheales+Polypodiales) to exclusion of Salviniales
tar_file(
  plastome_order_constraint_file,
  write_order_constraint_tree(
    plastome_alignment,
    ppgi_taxonomy,
    path(int_dir, "iqtree/plastome/order_constraint.tre")
  )
),
# Full plastome tree: all species, with 3rd positions, constrained by backbone
tar_target(
  plastome_tree,
  iqtree(
    plastome_alignment,
    spp = plastome_partition_file,
    m = "MFP",
    bb = 1000,
    nt = plastome_tree_nt_setting,
    seed = 20220123,
    redo = plastome_tree_redo_setting,
    echo = TRUE,
    wd = path(int_dir, "iqtree/plastome"),
    g = plastome_order_constraint_file,
    other_args = c("-t", "PARS"),
    tree_path = path(
      int_dir, "iqtree/plastome/plastome_partitions.txt.contree"
    )
  ),
  deployment = "main"
),
```

Even directly forcing the desired topology as a hard constraint did not
produce a satisfactory result (constraining to a relationship not otherwise
supported by the data is itself a strong signal that the underlying support
for the "expected" topology, at least from this plastome dataset alone, may
be weak or absent) — which is what ultimately motivated concluding this may
be a genuine feature of the data rather than a fixable artifact.

## Conclusion

None of the eleven approaches above — spanning partition scheme, model
selection (including FreeRate models), 3rd-codon-position removal via three
successive codon-aware alignment implementations, and three different
backbone-constraint strategies (self-derived exemplar tree, an independent
published tree, and a direct hard topological constraint) — resolved the
Cyatheales-sister-to-Salviniales placement in the plastome backbone tree.

Given how resistant this result was across such a wide range of
methodological changes, including a direct topological constraint, the
current working hypothesis is that this reflects genuine signal in FTOL's
plastome dataset rather than a correctable analysis artifact. This is
potentially of independent scientific interest and a candidate topic for a
focused follow-up paper examining this node specifically (e.g. site-specific
saturation analysis, gene-tree discordance, alternative topology tests). If
that investigation happens, this document and the commit history it
references (commits `b7e021e` through `05785d4`, plus the uncommitted
order-constraint code preserved above) are the starting point.

The codebase itself has been reverted to the pre-`b7e021e` non-partitioned
analysis method for both the plastome and Sanger trees so the standard FTOL
update pipeline could proceed with the next GenBank release.
