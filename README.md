# FTOL

[![DOI](https://zenodo.org/badge/475787005.svg)](https://zenodo.org/badge/latestdoi/475787005)

Code repository to generate an automatically updated **Fern Tree of Life**.

Please see the accompanying paper:
- Nitta JH, Schuettpelz E, Ramírez-Barahona S, Iwasaki W. An open and continuously updated fern tree of life (FTOL) https://doi.org/10.1101/2022.03.31.486640 

All code is in **R**, and workflow is controlled with the [targets](https://github.com/ropensci/targets) package.

## Docker

A [docker image](https://hub.docker.com/r/joelnitta/ftol) is available to run the code.

Docker tags match the version of FTOL; e.g., the image with tag `1.0.1` is the image used to generate FTOL v1.0.1.

## Setup

Run the [setup.R](R/setup.R) script to create the folder structure needed to store external files and download most of the data files automatically.

```
source("R/setup.R")
```

Alternatively, you can manually create the following folder hierarchy yourself:

```
_targets
└── user
    ├── data_raw
    │   ├── ref_aln
    │   └── restez
    │       └── sql_db
    ├── intermediates
    │   ├── blast_sanger
    │   ├── iqtree
    │   │   ├── plastome
    │   │   ├── sanger
    │   │   ├── sanger_1
    │   │   ├── sanger_10
    │   │   ├── sanger_2
    │   │   ├── sanger_3
    │   │   ├── sanger_4
    │   │   ├── sanger_5
    │   │   ├── sanger_6
    │   │   ├── sanger_7
    │   │   ├── sanger_8
    │   │   ├── sanger_9
    │   │   └── sanger_fast
    │   ├── ref_seqs
    │   └── treepl
    │       ├── con
    │       ├── ml
    │       └── ts
    └── results
```

Another folder called `ftol_data` (to store data files generated by this workflow that will be made available via the [ftolr R package](https://github.com/fernphy/ftolr)) also needs to be created in the project root. This folder is itself a repo that can be cloned from https://github.com/fernphy/ftol_data.

## Data

If `setup.R` was run successfully, it will have already downloaded and unzipped the input data files from FigShare.

Alternatively, you can do so manually following these instructions:

1. Download `ref_aln.tar.gz` (reference alignments) and `restez_sql_db.tar.gz` (local GenBank database) from FigShare (https://doi.org/10.6084/m9.figshare.19474316)
2. Unzip `ref_aln.tar.gz` and put the `ref_aln` folder in `_targets/user/data_raw/`
3. Unzip `restez_sql_db.tar.gz` and put the `bat` and `sql_logs` folders in `_targets/user/data_raw/restez/sql_db`.

## Running the code

### Data preparation scripts

- The script to generate the local GenBank database (`restez_sql_db.tar.gz`) is [setup_gb.R](R/setup_gb.R).

- The script (targets workflow) to generate the reference FASTA files for extracting target gene regions (`ref_aln.tar.gz`) is [prep_ref_seqs_plan.R](prep_ref_seqs_plan.R)

The data that result from these scripts have been made available on FigShare as described [above](#data), so these generally shouldn't need to be run.

### Main workflow

- [_targets.R](_targets.R) defines the main workflow to generate FTOL. This can be run with `targets::tar_make()`.

Note that this code was designed to be run on a multi-core machine, so the number of cores specified (e.g., [here](https://github.com/fernphy/ftol/blob/1c7569eb3bbd93864016bbc1b1df1d11f8d4d62c/_targets.R#L459)) may need to be changed.

The complete workflow takes 1-2 weeks to complete, with phylogenetic analysis taking up by far most of the time.

### Running with Docker

Launch a container in the background:

```
docker run \
  --rm \
  -dt \
  -v ${PWD}:/wd \
  -w /wd \
  -e USERID=$(id -u) \
  -e GROUPID=$(id -g) \
  joelnitta/ftol:latest bash
```

## License

- code: [MIT](LICENSE)
