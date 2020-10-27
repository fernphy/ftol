#!/bin/bash

# Example input argument
# arg="UFG_393201_P03_WA05:atpF"

# Split input into sample and gene, separated by colon
SAMPLE=$(echo $1 | cut -d ":" -f 1)
GENE=$(echo $1 | cut -d ":" -f 2)

# Extract name of reference bait matched by BLASTX
# (e.g., MH173078-psaA)
REF=`cat ~/ftol/intermediates/hybpiper/$SAMPLE/$GENE/${GENE}_baits.fasta | head -n 1 | grep -o '[^>]*'`

# Extract path to interleaved fasta file
# (path formatted for running inside docker container)
FASTA=/data/intermediates/hybpiper/$SAMPLE/$GENE/${GENE}_interleaved.fasta

# Format name of SAM output
SAM=${SAMPLE}-${GENE}.sam

# Format name of consenus FASTA output
CON=${SAMPLE}-${GENE}-con.fasta

# Run BBMAP, direct output to log named by sample + gene
docker run --rm -w /data/intermediates/bbmap --user root \
  -v /home/joelnitta/ftol/:/data \
  quay.io/biocontainers/bbmap:38.87--h1296035_0 bbmap.sh \
  in=$FASTA \
  out=$SAM \
  t=1 \
  ref=/data/intermediates/ref_dna/$REF.fasta \
  \nodisk 2> ./intermediates/bbmap/${SAMPLE}-${GENE}-log.txt

# activate kindel
source ~/miniconda2/etc/profile.d/conda.sh
conda activate kindel

# Extract consensus from SAM with min-depth of 1, trim Ns from alignment ends
kindel consensus ~/ftol/intermediates/bbmap/$SAM --min-depth 1 -t >  ~/ftol/intermediates/bbmap/$CON

# Remove temporary files
rm -f ~/ftol/intermediates/bbmap/$SAM
