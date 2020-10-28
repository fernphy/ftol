#!/bin/bash

# Run align_hybpip_bbmap.sh in parallel over 20 cores
# requires gnu parallel to be installed.
parallel -a working/align_hybpiper_args.txt \
  --eta \
  -j 20 \
  --max-args 1 \
  --load 80% \
  --noswap \
  R/align_hybpip_bbmap.sh
