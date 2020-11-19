# Raw Data README

## Raw data files

Raw data files used in this analysis.

- `goflag`: Folder containing raw reads of fern pilot samples for the
  [GoFlag](http://flagellateplants.group.ufl.edu/) project (48 samples total).
  Downloaded from [Globus](https://www.globus.org/) file transfer program on
  2020-08-20. For a list of MD5 hashes for each raw data file (paired-end
  fastq.gz files), see `md5_checksums.csv`. Includes one set of reads after target
  enrichment (`pilot_raw_reads`), and one set of reads without enrichment
  (`skimming_raw_reads`). For details on enrichment methods, see [GoFlag pilot
  pre-print](https://doi.org/10.1101/2020.05.29.124081).

- `genbank_names_with_mult_syns_select.csv`: List of manually selected names for
  taxa with Sanger sequences downloaded from GenBank that had multiple synonyms
  detected during taxonomic name resolution that could not be determined
  automatically. Column `query` indicates the queried (original) name;
  `name_resolved_manual` is the manually selected name to use for analysis.
  Comma-separated text file.

- `Pilot_Ferns_TargetCapture_Skimming.txt`: Metadata of fern [pilot
  samples](https://doi.org/10.1101/2020.05.29.124081) for the
  [GoFlag](http://flagellateplants.group.ufl.edu/) project. Columns: `Taxon`
  (taxon name), `Herbarium` (acronym of herbarium where voucher specimen is
  lodged), `Voucher` (voucher specimen), `NCBI Accession - Targeted Capture`
  (NCBI accession number for target capture data), `Targeted Capture ID` (ID for
  target capture reads), `Genome Skimming ID` (ID for genome skimming reads).
  Tab-separated text file.

- `plastid_targets.RDS`: Temporary R data file with plastid gene sequences used as
  targets for [HybPiper](https://github.com/mossmatters/HybPiper). This will be
  generated as part of the main analysis workflow.

- `plastome_outgroups.csv`: Plastome accessions used as outgroups in plastid
  phylogenetic analysis. Columns: `group name` (major plant group), `species`
  (species name), `acccession` (GenBank accession number). Comma-separated text 
  file.

- `ppgi_taxonomy_mod.csv`: Taxonomic system for ferns and lycophytes at genus
  level and higher following Pteridophyte Phylogeny Group (PPG) I 2016, modified
  with new genera added since that time and additional genera used in the World
  Ferns list but not included in PPG I. Comma-separated text file.

- `taxa.txt`: Taxonomic data of tracheophytes from the Species 2000 & ITIS
  Catalogue of Life, 2019 edition. Tab-separated file including 1,168,025
  observations of 31 variables, mostly taxonomic ranks and names in the Darwin
  Core Archive format. Originally included in a zip file downloaded from
  http://www.catalogueoflife.org/DCA_Export/zip/archive-kingdom-plantae-phylum-tracheophyta-bl3.zip
  on 2019-05-19. Data have been subset to tracheophytes by selecting "Plantae"
  for "Top level group" and "Tracheophyta" for "Phylum" at
  http://www.catalogueoflife.org/DCA_Export/index.php. Names of pteridophytes in
  this file come from the [World Ferns
  list](https://worldplants.webarchiv.kit.edu/ferns/) Nov 2018 version by
  Michael Hassler.

- `testo_sundue_2016_calibrations.csv`: Calibration dates for dating fern
  phylogeny based on Testo and Sundue (2016). Columns: `Clade` (name of clade),
  `Stem_Crown` (status of clade as stem group or crown group), `Fossil` (name of
  fossil), `Age` (fossil age in millions of years before present), `Age_type`
  (age type used by treepl), `Citation` (reference for the fossil), `taxon_1`
  and `taxon_2` are names of tips in the tree that have a most recent common
  ancestor corresponding to that clade. Comma-separated text file.

- `wei_2017_coding_genes.txt`: List of coding genes (83 gene set) used in Wei et 
  al. (2017).

## Additional documentation

These files are not used in analysis, but provide additional documentation.

- `md5_checksums.csv`: List of MD5 checksums for each raw data file. Columns:
  `target` (name of target corresponding to the file path in the drake plan),
  `file` (file name), `path` (path to the file), `copy` (TRUE/FALSE indicating
  whether the file should be copied when creating versioned data archives),
  `hash` (MD5 checksum). Only raw data files that are subject to change are
  copied into versioned data archives. Comma-separated text file.

## References

Testo, Weston, and Michael Sundue. "A 4000-species dataset provides new insight
into the evolution of ferns." Molecular Phylogenetics and Evolution 105 (2016):
200-211

Wei, Ran, et al. "Plastid phylogenomics resolve deep relationships among
eupolypod II ferns with rapid radiation and rate heterogeneity." Genome biology
and evolution 9.6 (2017): 1646-1657.
