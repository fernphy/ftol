This README.txt file was generated on 2023-08-15 by Joel Nitta

--------------------------------------------------------------------------------

GENERAL INFORMATION

--------------------------------------------------------------------------------

Title of Dataset: Fern Tree of Life (FTOL) input data

Principal Investigator: Joel H. Nitta

Department of Integrated Biosciences, Graduate School of Frontier Sciences, The
University of Tokyo, Chiba, Japan. joelnitta@gmail.com

Associate or Co-investigators: Eric Schuettpelz, Santiago Ramírez-Barahona,
Wataru Iwasaki

Date of data collection: 1990 - 2022

Geographic location of data collection: Global

Information about funding sources that supported the collection of the data:
Funding provided in part by the Japan Society for the Promotion of Science
(Kakenhi) Grant Number 16H06279 and the Smithsonian National Museum of Natural
History Peter Buck Fellowship (JHN).

--------------------------------------------------------------------------------

SHARING/ACCESS INFORMATION

--------------------------------------------------------------------------------

Licenses/restrictions placed on the data: CC0 v1.0 license

Links to publications that cite or use the data:

Nitta JH, Schuettpelz E, Ramírez-Barahona S, Iwasaki W. 2022. An open and
continuously updated fern tree of life (FTOL).
https://doi.org/10.1101/2022.03.31.486640

Links to other publicly accessible locations of the data: none.

Links/relationships to ancillary data sets:

-   ferncal (https://github.com/fernphy/ferncal)
-   pteridocat (https://github.com/fernphy/pteridocat)
-   FTOL input data on FigShare
    (https://doi.org/10.6084/m9.figshare.19474316.v1)

Was data derived from another source?

-   ppgi_taxonomy_mod.csv is derived from PPGI.csv
    (https://gist.github.com/mutolisp/d0bd8e7c998303e316ae09e4f37c91ae)
-   equisetum_subgenera.csv is based on data from
    https://en.wikipedia.org/wiki/Equisetum
-   target_coding_genes.txt is based on data in Wei et al. (2017)

Recommended citation for this dataset:

Nitta JH, Schuettpelz E, Ramírez-Barahona S, Iwasaki W. 2022. An open and
continuously updated fern tree of life (FTOL).
https://doi.org/10.1101/2022.03.31.486640

--------------------------------------------------------------------------------

DATA & FILE OVERVIEW

--------------------------------------------------------------------------------

File List:

-   accs_exclude.csv: GenBank accessions to exclude from analysis.
-   accs_manual.csv: GenBank accessions to include in analysis.
-   equisetum_subgenera.csv: Equisteum subgenera.
-   plastome_outgroups.csv: GenBank accessions of non-fern plastome sequences to
    include in analysis.
-   ppgi_taxonomy_mod.csv: Modified Pteridophyte Phylogeny Group I taxonomic
    system.
-   target_coding_genes.txt: List of coding genes.

--------------------------------------------------------------------------------

METHODOLOGICAL INFORMATION

--------------------------------------------------------------------------------

The data included here are used in a pipeline that (mostly) automatically
generates a maximally sampled fern phylogenetic tree based on plastid sequences
in GenBank (https://github.com/fernphy/ftol).

The first step is to generate a set of reference FASTA files for 79 target loci
(one per locus; ref_aln.tar.gz). These include 77 protein-coding genes based on
a list of 83 genes (Wei et al. 2017) that was filtered to only genes that show
no evidence of duplication (target_coding_genes.txt), plus two spacer regions
(trnL-trnF and rps4-trnS). Each FASTA file in ref_aln.tar.gz includes one
representative (longest) sequence per avaialable fern genus. This is done with
custom R scripts contained in https://github.com/fernphy/ftol, in particular
prep_ref_seqs_plan.R
(https://github.com/fernphy/ftol/blob/main/prep_ref_seqs_plan.R).

Next, all available fern accessions for seven target “Sanger loci” (plastid
regions typically sequenced using Sanger technology) and all available fern
plastomes (accessions >7000 bp) are downloaded from GenBank. Non-fern accessions
listed in plastome_outgroups.csv are downloaded as well. Sequences matching the
target loci are then extracted from each accesion using the FASTA files
contained in ref_aln.tar.gz as references with the “Reference_Blast_Extract.py”
script of superCRUNCH (Portik and Wiens 2020). Any accessions matching those
listed in accs_exclude.csv are excluded as putative rogues (i.e.,
misidentifications or contaminations). Any accessions matching those listed in
accs_include.csv are used regardless of other accessions for the same species in
GenBank.

The extracted sequences are aligned with MAFFT (Katoh et al. 2002), phylogenetic
analysis is done using IQ-TREE (Nguyen et al. 2015) and divergence times
estimated with treePL (Smith and O’Meara 2012). During molecular dating,
equisetum_subgenera.csv is used to specify some clades within Equisetum whose
ages are constrained by fossils, and ppgi_taxonomy_mod.csv is used to map
higher-level clade names (e.g., family, order, etc.) to species (tips of the
phylogeny).

For additional methodological details, see:

Nitta JH, Schuettpelz E, Ramírez-Barahona S, Iwasaki W. 2022. An open and
continuously updated fern tree of life (FTOL).
https://doi.org/10.1101/2022.03.31.486640

--------------------------------------------------------------------------------

DATA-SPECIFIC INFORMATION

--------------------------------------------------------------------------------

accs_exclude.csv: GenBank accessions to exclude from analysis.

Number of variables: 3

Number of cases/rows: 58

Variable list:

-   accession: GenBank accession number.
-   gene: Name of locus or if whole plastome, ‘plastome’.
-   note: Note about why the accession is excluded.

Missing data codes: No missing data.

Specialized formats or other abbreviations used: None.

MD5 checksum: ee64b3cd3dc22ed1d25b00e08e9d1d19

--------------------------------------------------------------------------------

accs_manual.csv: GenBank accessions to preferentially include in analysis.
Accessions in this list will be used regardless of other accessions for the same
species in GenBank.

Number of variables: 9

Number of cases/rows: 4

Variable list:

-   species: Species name (same as in FTOL tree).
-   atpA: atpA GenBank accession.
-   atpB: atpB GenBank accession.
-   matK: matK GenBank accession.
-   rbcL: rbcL GenBank accession.
-   rps4: rps4 GenBank accession.
-   rps4-trnS: rps4-trnS GenBank accession.
-   trnL-trnF: trnL-trnF GenBank accession.
-   note: Note on usage

Missing data codes: Non-applicable data indicated with ‘NA’.

Specialized formats or other abbreviations used: None.

MD5 checksum: 472dffe6f6bdcd5069f83d351c705e1a

--------------------------------------------------------------------------------

equisetum_subgenera.csv: Equisteum subgenera

Number of variables: 2

Number of cases/rows: 47

Variable list:

-   scientificName: Scientific name.
-   subgenus: Subgenus.

Missing data codes: Nothing entered if there is no assigned subgenus.

Specialized formats or other abbreviations used: None.

MD5 checksum: 1b964b8ff9de5557d52a7f3047f7ab5a

--------------------------------------------------------------------------------

plastome_outgroups.csv: GenBank accessions of non-fern plastome sequences to
include in analysis.

Number of variables: 3

Number of cases/rows: 19

Variable list:

-   group: Name of clade containing the species.
-   species: Species.
-   accession: GenBank accession number.

Missing data codes: None.

Specialized formats or other abbreviations used: None.

MD5 checksum: 93f08195e9dbd27d0b84a19d05180729

--------------------------------------------------------------------------------

ppgi_taxonomy_mod.csv: Taxonomic classification of ferns and lycophytes based on
Pteridophyte Phylogeny Group I (2016), with modifications.

Number of variables: 8

Number of cases/rows: 407

Variable list:

-   class: Class.
-   order: Order.
-   suborder: Suborder.
-   family: Family.
-   subfamily: Subfamily.
-   genus: Genus.
-   notes: Notes, in particular if the taxon differs from PPG I (2016).
-   nothogenus: Is the taxon a nothogenus? ‘yes’ or ‘no’.

Missing data codes: Non-applicable data indicated with ‘NA’.

Specialized formats or other abbreviations used: None.

MD5 checksum: d1612ac49c439b21e848e542db706e72

--------------------------------------------------------------------------------

CHANGE LOG

2022-05-27

-   Move csv and plain text files from FigShare to Github.

2022-04-04

-   Update DOIs.

2022-03-31

-   Generate this README file.

--------------------------------------------------------------------------------

REFERENCES

Katoh, Kazutaka, Kazuharu Misawa, Keiichi Kuma, and Takashi Miyata. 2002.
“MAFFT: A Novel Method for Rapid Multiple Sequence Alignment Based on Fast
Fourier Transform.” Nucleic Acids Research 30 (14): 3059–66.
https://doi.org/10.1093/nar/gkf436.

Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh.
2015. “IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating
Maximum-Likelihood Phylogenies.” Molecular Biology and Evolution 32 (1): 268–74.
https://doi.org/10.1093/molbev/msu300.

Portik, Daniel M., and John J. Wiens. 2020. “SuperCRUNCH: A Bioinformatics
Toolkit for Creating and Manipulating Supermatrices and Other Large Phylogenetic
Datasets.” Edited by David Orme. Methods in Ecology and Evolution 11 (6):
763–72. https://doi.org/ggx588.

Pteridophyte Phylogeny Group I. 2016. “A Community-Derived Classification for
Extant Lycophytes and Ferns.” Journal of Systematics and Evolution 54 (6):
563–603. https://doi.org/10.1111/jse.12229.

Smith, Stephen A, and Brian C. O’Meara. 2012. “treePL: Divergence Time
Estimation Using Penalized Likelihood for Large Phylogenies.” Bioinformatics 28
(20): 2689–90. https://doi.org/10.1093/bioinformatics/bts492.

Wei, Ran, Yue-Hong Yan, AJ Harris, Jong-Soo Kang, Hui Shen, Qiao-Ping Xiang, and
Xian-Chun Zhang. 2017. “Plastid Phylogenomics Resolve Deep Relationships Among
Eupolypod II Ferns with Rapid Radiation and Rate Heterogeneity.” Genome Biology
and Evolution 9 (6): 1646–57. https://doi.org/10.1093/gbe/evx107.

Katoh, Kazutaka, Kazuharu Misawa, Keiichi Kuma, and Takashi Miyata. 2002.
“MAFFT: A Novel Method for Rapid Multiple Sequence Alignment Based on Fast
Fourier Transform.” Nucleic Acids Research 30 (14): 3059–66.
https://doi.org/10.1093/nar/gkf436.

Nguyen, Lam-Tung, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh.
2015. “IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating
Maximum-Likelihood Phylogenies.” Molecular Biology and Evolution 32 (1): 268–74.
https://doi.org/10.1093/molbev/msu300.

Portik, Daniel M., and John J. Wiens. 2020. “SuperCRUNCH: A Bioinformatics
Toolkit for Creating and Manipulating Supermatrices and Other Large Phylogenetic
Datasets.” Edited by David Orme. Methods in Ecology and Evolution 11 (6):
763–72. https://doi.org/ggx588.

Pteridophyte Phylogeny Group I. 2016. “A Community-Derived Classification for
Extant Lycophytes and Ferns.” Journal of Systematics and Evolution 54 (6):
563–603. https://doi.org/10.1111/jse.12229.

Smith, Stephen A, and Brian C. O’Meara. 2012. “treePL: Divergence Time
Estimation Using Penalized Likelihood for Large Phylogenies.” Bioinformatics 28
(20): 2689–90. https://doi.org/10.1093/bioinformatics/bts492.

Wei, Ran, Yue-Hong Yan, AJ Harris, Jong-Soo Kang, Hui Shen, Qiao-Ping Xiang, and
Xian-Chun Zhang. 2017. “Plastid Phylogenomics Resolve Deep Relationships Among
Eupolypod II Ferns with Rapid Radiation and Rate Heterogeneity.” Genome Biology
and Evolution 9 (6): 1646–57. https://doi.org/10.1093/gbe/evx107.
