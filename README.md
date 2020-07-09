## FTOL

The goal of this project is to generate an automatically updated **Fern Tree of Life**.

All code is in **R**, and workflow is controlled with the [drake](https://docs.ropensci.org/drake/) package.

### Docker

A [docker image](https://hub.docker.com/repository/docker/joelnitta/ftol) is available to run the code.

### Data

Raw data files should be downloaded from Dropbox using the following links, and saved to `data_raw` in the project root:

- [Catalog of Life taxonomic data for tracheophytes (`archive-kingdom-plantae-phylum-tracheophyta-bl3/taxa.txt`)](https://www.dropbox.com/s/uk49g48h8jhoslv/taxa.txt?dl=0)

- [Modified PPGI taxonomy (`ppgi_taxonomy_mod.csv`)](https://www.dropbox.com/s/l5ptrndae4jy53q/ppgi_taxonomy_mod.csv?dl=0)

- [List of plastid coding genes from Wei et al 2017 (`wei_2017_coding_genes.txt`)](https://www.dropbox.com/s/aozhk47lguzfxnj/wei_2017_coding_genes.txt?dl=0)

- [Outgroup plastome accessions (`plastome_outgroups.csv`)](https://www.dropbox.com/s/b2m1ln9e5itjw43/plastome_outgroups.csv?dl=0)

- [Calibration dates after Testo and Sundue 2016 (`testo_sundue_2016_calibrations.csv`)](https://www.dropbox.com/s/mg1k8zzlwko24bf/testo_sundue_2016_calibrations.csv?dl=0)

- [Manually selected synonyms for resolving names of plastid genes (`genbank_names_with_mult_syns_select.csv`)](https://www.dropbox.com/s/8qoae4kl8ps1z6h/genbank_names_with_mult_syns_select.csv?dl=0)
 
### Running the code

To run the analysis, execute `drake::r_make()` from the root of the repo.

### Details

For more info, [see the wiki](https://github.com/joelnitta/ftol/wiki).
