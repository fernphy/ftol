#!/bin/bash

# This aligns short reads that were sorted by HybPiper with BLASTX against
# reference "target" gene sequences. Assumes hybpiper has been run with output
# to $WDIR/intermediates/hybpiper, and all reference sequences available
# (one gene per fasta file) in $WDIR/intermediates/ref_dna.

# Will output results to $WDIR/intermediates/bbmap

# Takes a single input argument.
# Example: (sample name, gene name, separated by colon)
# UFG_393201_P03_WA05:atpF

# Set and enter working directory
WDIR=/home/joelnitta/ftol
cd $WDIR

# Split input into sample and gene, separated by colon
SAMPLE=$(echo $1 | cut -d ":" -f 1)
GENE=$(echo $1 | cut -d ":" -f 2)

# Extract name of reference bait matched by BLASTX in hybpiper
# (e.g., MH173078-psaA)
REF=`cat $WDIR/intermediates/hybpiper/$SAMPLE/$GENE/${GENE}_baits.fasta | head -n 1 | grep -o '[^>]*'`

# Extract path to interleaved fasta file
# (path formatted for running inside docker container)
FASTA=/data/hybpiper_sample/$GENE/${GENE}_interleaved.fasta

# Format name of SAM output
SAM=${SAMPLE}-${GENE}.sam

# Format name of consenus FASTA output
CON=${SAMPLE}-${GENE}.fasta

# This will create a lot of files (two per sample/gene combination) when we loop.
# To avoid hitting docker's file number limit, only mount subfolders.
# (don't issue an error, and don't create a new folder, if one already exists)
mkdir -p $WDIR/intermediates/bbmap/${SAMPLE}-${GENE}

# Map reads to reference with BBMAP, direct output to log named by sample + gene.
# Outputs in SAM format.
# -t: number of threads
# -Xmx5g: use 1g memory per job (java argument)
# \nodisk: don't write the index files to disk
# bbmap prints to stderr, so to log output use `2>`
docker run --rm -w /data/bbmap/${SAMPLE}-${GENE} --user root \
  -v $WDIR/intermediates/hybpiper/$SAMPLE:/data/hybpiper_sample \
  -v $WDIR/intermediates/ref_dna:/data/ref_dna \
  -v $WDIR/intermediates/bbmap/${SAMPLE}-${GENE}:/data/bbmap/${SAMPLE}-${GENE} \
  quay.io/biocontainers/bbmap:38.87--h1296035_0 bbmap.sh \
  in=$FASTA \
  out=$SAM \
  t=1 \
  -Xmx1g \
  ref=/data/ref_dna/$REF.fasta \
  \nodisk 2> $WDIR/intermediates/bbmap/${SAMPLE}-${GENE}/bbmap.log

# Activate kindel
# assumes miniconda2 has been installed, and a conda environment
# called 'kindel' already created
source ~/miniconda2/etc/profile.d/conda.sh
conda activate kindel

# Extract consensus from SAM with min-depth of 1, trim Ns from alignment ends
kindel consensus $WDIR/intermediates/bbmap/${SAMPLE}-${GENE}/$SAM --min-depth 1 -t \
  > $WDIR/intermediates/bbmap/${SAMPLE}-${GENE}/$CON \
  2> $WDIR/intermediates/bbmap/${SAMPLE}-${GENE}/kindel.log
