#!/bin/bash
# run this once per day at midnight by adding to crontab (with crontab -e):
# 0 0 * * * bash ~/ftol/setup_gb.sh > ~/ftol/setup_gb.log
cd ~/ftol
docker run --rm -v ${PWD}:/wd -w /wd \
  joelnitta/ftol:latest Rscript -e "R/setup_gb.R"
