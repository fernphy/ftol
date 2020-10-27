#!/bin/bash

N=10
(
while read args; do
  ((i=i%N)); ((i++==0)) && wait
  R/align_hybpip_bbmap.sh "$args" &
done < working/align_hybpiper_args.txt
)
