# FTOL

Code repository to generate an automatically updated **Fern Tree of Life**.

Please see the accompanying paper: Nitta JH, Schuettpelz E, Ramírez-Barahona S, Iwasaki W. An open and continuously updated fern tree of life (FTOL)

All code is in **R**, and workflow is controlled with the [targets](https://github.com/ropensci/targets) package.

## Docker

A [docker image](https://hub.docker.com/r/joelnitta/ftol) is available to run the code.

## Setup

First run `targets::tar_make()` from the root of the repo. This won't finish (since there's no data yet!) but it will set up the folders needed for storing various files.

You should now have a folder called `_targets` present in the repo.

Next, run the [setup.R](R/setup.R) script to create all the subfolders needed to store external files in `_targets`. It also downloads most of the needed data files automatically.

```
source("R/setup.R")
```

## Data

### FTOL input data on Figshare

If `setup.R` was run successfully, it will have already downloaded and unzipped the needed files from FigShare.

Alternatively, you can do so manually following these instructions:

Download all of the files contained in the [FTOL input data on FigShare](https://doi.org/10.6084/m9.figshare.19474316.v1) into `_targets/user/data_raw` by clicking on the "Download All" button.

Unzip the files and place them in the `_targets/user/data_raw` folder. You will also need to unzip the `ref_aln.tar.gz` archive that is contained within the zip folder.

### Other data files

Two other data files need to be manually downloaded and placed in `_targets/user/data_raw`:

- Supplemental Data 1  (1-s2.0-S1055790316302287-mmc2.xlsx) from [Testo and Sundue (2016) Mol. Phylogenetics Evol.](https://doi.org/10.1016/j.ympev.2016.09.003)
- NCBI taxonomy database dump 2022-02-01 (taxdmp_2022-02-01.zip). This can be downloaded from the [NCBI FTP server](https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2022-02-01.zip)
## Running the code

To run the analysis, run `targets::tar_make()` from the root of the repo.