#!/bin/bash

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

# Make bowtie2 index file from a single fasta reference
docker run --rm -w /data/intermediates/bowtie2 --user root \
  -v /home/joelnitta/ftol/:/data \
  biocontainers/bowtie2:v2.4.1_cv1 bowtie2-build \
  /data/intermediates/ref_dna/$REF.fasta \
  $REF
  
# Align reads to reference using bowtie2
# -x: name of reference index
# -f: flag that input is in fasta format
# --interleaved: flag that input is interleaved
# -S: name of output SAM file
docker run --rm -w /data/intermediates/bowtie2 --user root \
  -v /home/joelnitta/ftol/:/data \
  biocontainers/bowtie2:v2.4.1_cv1 bowtie2 \
  -x $REF \
  -f \
  --interleaved $FASTA \
  --very-sensitive \
  -S $SAM
  
# activate kindel
source ~/miniconda2/etc/profile.d/conda.sh
conda activate kindel

# Extract consensus from SAM with min-depth of 1
kindel consensus ~/ftol/intermediates/bowtie2/$SAM --min-depth 1 >  ~/ftol/intermediates/bowtie2/$CON

# Remove temporary files
rm ~/ftol/intermediates/bowtie2/$REF.*
rm ~/ftol/intermediates/bowtie2/$SAM

