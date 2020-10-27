#!/bin/bash

cat working/align_hybpiper_args.txt | xargs -P10 -n1 bash -c 'R/align_hybpip_bbmap.sh "$1"' {}
