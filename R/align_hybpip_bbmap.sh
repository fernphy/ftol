#!/bin/bash

# This aligns short reads that were sorted by HybPiper with BLASTX against
# reference "target" gene sequences. Assumes hybpiper has been run with output
# to ./intermediates/hybpiper, and all reference sequences available
# (one gene per fasta file) in ./intermediates/ref_dna.

# Example input argument (sample name, gene name, separated by colon)
# arg="UFG_393201_P03_WA05:atpF"

# Working directory is /home/joelnitta/ftol/

# Split input into sample and gene, separated by colon
SAMPLE=$(echo $1 | cut -d ":" -f 1)
GENE=$(echo $1 | cut -d ":" -f 2)

# Extract name of reference bait matched by BLASTX in hybpiper
# (e.g., MH173078-psaA)
REF=`cat ~/ftol/intermediates/hybpiper/$SAMPLE/$GENE/${GENE}_baits.fasta | head -n 1 | grep -o '[^>]*'`

# Extract path to interleaved fasta file
# (path formatted for running inside docker container)
FASTA=/data/intermediates/hybpiper/$SAMPLE/$GENE/${GENE}_interleaved.fasta

# Format name of SAM output
SAM=${SAMPLE}-${GENE}.sam

# Format name of consenus FASTA output
CON=${SAMPLE}-${GENE}-con.fasta

# Map reads to reference with BBMAP, direct output to log named by sample + gene.
# Outputs in SAM format.
# -t: number of threads
# -Xmx5g: use 5g memory per job (java argument)
# \nodisk: don't write the index files to disk
# bbmap prints to stderr, so to log output use `2>`
docker run --rm -w /data/intermediates/bbmap --user root \
  -v /home/joelnitta/ftol/:/data \
  quay.io/biocontainers/bbmap:38.87--h1296035_0 bbmap.sh \
  in=$FASTA \
  out=$SAM \
  t=1 \
  -Xmx5g \
  ref=/data/intermediates/ref_dna/$REF.fasta \
  \nodisk 2> ./intermediates/bbmap/${SAMPLE}-${GENE}-log.txt

# Activate kindel
source ~/miniconda2/etc/profile.d/conda.sh
conda activate kindel

# Extract consensus from SAM with min-depth of 1, trim Ns from alignment ends
kindel consensus ~/ftol/intermediates/bbmap/$SAM --min-depth 1 -t >  ~/ftol/intermediates/bbmap/$CON

# Remove SAM file, since we only want the consensus FASTA
rm -f ~/ftol/intermediates/bbmap/$SAM
